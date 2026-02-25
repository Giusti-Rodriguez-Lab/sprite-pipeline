use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use clap::Parser;
use crossbeam_channel::bounded;
use rayon::prelude::*;

use barcode_id::config::Config;
use barcode_id::tag::TagMaps;
use barcode_id::{fastq, matcher};

const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Records per batch handed to the rayon thread pool.
const BATCH_SIZE: usize = 100_000;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

/// Fast SPRITE barcode identification.
///
/// Drop-in replacement for BarcodeIdentification_v1.2.0.jar.
/// Reads paired FASTQ.gz files, identifies barcodes in each read pair using
/// the provided config, and writes barcoded FASTQ.gz files where every read
/// name has the form:  @readname::[Tag1][Tag2]...[NOT_FOUND]...
#[derive(Parser, Debug)]
#[command(name = "barcode_id", version = VERSION)]
struct Args {
    /// Input FASTQ R1 (gzipped)
    #[arg(long)]
    input1: PathBuf,

    /// Input FASTQ R2 (gzipped)
    #[arg(long)]
    input2: PathBuf,

    /// Output FASTQ R1 (gzipped)
    #[arg(long)]
    output1: PathBuf,

    /// Output FASTQ R2 (gzipped)
    #[arg(long)]
    output2: PathBuf,

    /// Barcode configuration file
    #[arg(long)]
    config: PathBuf,

    /// Number of worker threads for parallel processing (0 = use all available CPUs)
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or("info"),
    )
    .init();

    let start = Instant::now();
    let args  = Args::parse();

    // --- Thread pool ---
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap_or_else(|e| {
                log::warn!("Could not configure thread pool size: {}", e);
            });
    }
    log::info!(
        "Using {} rayon threads.",
        rayon::current_num_threads(),
    );

    // --- Config and tag maps (Arc so rayon closures can borrow across threads) ---
    log::info!("Parsing configuration file: {:?}", args.config);
    let config = Arc::new(Config::from_file(&args.config).unwrap_or_else(|e| {
        log::error!("Failed to parse config: {}", e);
        std::process::exit(1);
    }));
    log::info!(
        "Loaded {} tag definitions. spacer_len={} laxity={}",
        config.tag_defs.len(),
        config.spacer_len,
        config.laxity,
    );

    log::info!("Building tag lookup maps...");
    let tag_maps = Arc::new(TagMaps::build(&config));
    log::info!("Tag maps ready.");

    // --- Channels (bound=2 for backpressure) ---
    // channel A: reader thread  → coordinator (raw paired batches)
    // channel B: coordinator    → writer thread (processed results)
    let (batch_tx, batch_rx)   = bounded::<Vec<(fastq::Record, fastq::Record)>>(2);
    let (result_tx, result_rx) = bounded::<Vec<(fastq::Record, fastq::Record)>>(2);

    // --- Reader thread: decompress + batch FASTQ pairs ---
    let input1 = args.input1;
    let input2 = args.input2;
    let reader_handle = std::thread::spawn(move || {
        let mut reader1 = fastq::open_gz(&input1).unwrap_or_else(|e| {
            log::error!("Cannot open input1 {:?}: {}", input1, e);
            std::process::exit(1);
        });
        let mut reader2 = fastq::open_gz(&input2).unwrap_or_else(|e| {
            log::error!("Cannot open input2 {:?}: {}", input2, e);
            std::process::exit(1);
        });

        loop {
            let mut batch: Vec<(fastq::Record, fastq::Record)> =
                Vec::with_capacity(BATCH_SIZE);
            let mut eof = false;

            for _ in 0..BATCH_SIZE {
                let r1 = match reader1.next_record() {
                    Ok(Some(r)) => r,
                    Ok(None) => { eof = true; break; }
                    Err(e) => {
                        log::error!("Read error on input1: {}", e);
                        std::process::exit(1);
                    }
                };
                let r2 = match reader2.next_record() {
                    Ok(Some(r)) => r,
                    Ok(None) => {
                        log::warn!("input2 ended before input1; output may be incomplete.");
                        eof = true;
                        break;
                    }
                    Err(e) => {
                        log::error!("Read error on input2: {}", e);
                        std::process::exit(1);
                    }
                };
                batch.push((r1, r2));
            }

            if batch.is_empty() {
                break;
            }
            if batch_tx.send(batch).is_err() {
                break; // coordinator dropped (shouldn't happen)
            }
            if eof {
                break;
            }
        }
        // batch_tx dropped here → coordinator's for-loop terminates
    });

    // --- Writer thread: write processed pairs to gzip output ---
    let output1 = args.output1;
    let output2 = args.output2;
    let writer_handle = std::thread::spawn(move || {
        let mut writer1 = fastq::create_gz(&output1).unwrap_or_else(|e| {
            log::error!("Cannot create output1 {:?}: {}", output1, e);
            std::process::exit(1);
        });
        let mut writer2 = fastq::create_gz(&output2).unwrap_or_else(|e| {
            log::error!("Cannot create output2 {:?}: {}", output2, e);
            std::process::exit(1);
        });

        let mut count: u64 = 0;
        for results in result_rx {
            for (out1, out2) in results {
                out1.write_to(&mut writer1).unwrap_or_else(|e| {
                    log::error!("Write error on output1: {}", e);
                    std::process::exit(1);
                });
                out2.write_to(&mut writer2).unwrap_or_else(|e| {
                    log::error!("Write error on output2: {}", e);
                    std::process::exit(1);
                });
                count += 1;
                if count % 500_000 == 0 {
                    log::info!("Processing read {}.", count);
                }
            }
        }

        // Flush gzip trailers — must call finish() or output files will be truncated.
        writer1.into_inner().unwrap().finish().unwrap();
        writer2.into_inner().unwrap().finish().unwrap();
        count
    });

    // --- Coordinator (main thread): receive batches → rayon → send results ---
    //
    // Timeline:
    //   reader:  [read][read][read]...
    //   rayon:         [proc][proc]...
    //   writer:              [writ][writ]...
    //
    // par_iter().collect() preserves batch-internal order; batch order is
    // preserved by the serial coordinator loop.
    for batch in batch_rx {
        let results: Vec<(fastq::Record, fastq::Record)> = batch
            .par_iter()
            .map(|(r1, r2)| matcher::process_pair(r1, r2, &config, &tag_maps))
            .collect();
        if result_tx.send(results).is_err() {
            log::error!("Writer thread disconnected unexpectedly.");
            std::process::exit(1);
        }
    }

    // Signal the writer that no more results are coming.
    drop(result_tx);

    reader_handle.join().expect("Reader thread panicked");
    let count = writer_handle.join().expect("Writer thread panicked");

    log::info!("Program complete.");
    log::info!("{} milliseconds elapsed.", start.elapsed().as_millis());
    log::info!("{} read pairs processed.", count);
}
