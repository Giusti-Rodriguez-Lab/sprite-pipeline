use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
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

    // --- Config ---
    log::info!("Parsing configuration file: {:?}", args.config);
    let config = Config::from_file(&args.config).unwrap_or_else(|e| {
        log::error!("Failed to parse config: {}", e);
        std::process::exit(1);
    });
    log::info!(
        "Loaded {} tag definitions. spacer_len={} laxity={}",
        config.tag_defs.len(),
        config.spacer_len,
        config.laxity,
    );

    // --- Tag maps ---
    log::info!("Building tag lookup maps...");
    let tag_maps = TagMaps::build(&config);
    log::info!("Tag maps ready.");

    // --- I/O ---
    let mut reader1 = fastq::open_gz(&args.input1).unwrap_or_else(|e| {
        log::error!("Cannot open input1 {:?}: {}", args.input1, e);
        std::process::exit(1);
    });
    let mut reader2 = fastq::open_gz(&args.input2).unwrap_or_else(|e| {
        log::error!("Cannot open input2 {:?}: {}", args.input2, e);
        std::process::exit(1);
    });
    let mut writer1 = fastq::create_gz(&args.output1).unwrap_or_else(|e| {
        log::error!("Cannot create output1 {:?}: {}", args.output1, e);
        std::process::exit(1);
    });
    let mut writer2 = fastq::create_gz(&args.output2).unwrap_or_else(|e| {
        log::error!("Cannot create output2 {:?}: {}", args.output2, e);
        std::process::exit(1);
    });

    // --- Batched parallel processing loop ---
    //
    // Each iteration:
    //   1. Read up to BATCH_SIZE pairs from the sequential gz readers (main thread).
    //   2. Process the batch in parallel via rayon — order is preserved by
    //      `par_iter` + `collect`, which reassembles results in input order.
    //   3. Write the ordered results to the output gz files (main thread).
    let mut count: u64 = 0;
    loop {
        // --- Read batch ---
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

        // --- Process batch in parallel (order-preserving) ---
        let results: Vec<(fastq::Record, fastq::Record)> = batch
            .par_iter()
            .map(|(r1, r2)| matcher::process_pair(r1, r2, &config, &tag_maps))
            .collect();

        // --- Write in order ---
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

        if eof {
            break;
        }
    }

    log::info!("Program complete.");
    log::info!("{} milliseconds elapsed.", start.elapsed().as_millis());
    log::info!("{} read pairs processed.", count);
}
