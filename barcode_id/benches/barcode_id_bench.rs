/// Criterion benchmarks for barcode_id.
///
/// Three benchmarks:
///   A. generate_neighbors  — Hamming expansion of a 15-mer at k=2
///   B. check_read          — full scSPRITE2 R2 layout on a realistic sequence
///   C. batch_throughput    — parallel processing of 1 000 pairs via rayon
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rayon::prelude::*;

use barcode_id::config::Config;
use barcode_id::fastq::Record;
use barcode_id::matcher::{check_read, process_pair};
use barcode_id::tag::{generate_neighbors, TagMaps};

// Production config embedded at compile time — no runtime file dependency.
const PROD_CONFIG: &str =
    include_str!("../../misc/config_dpm6_y-stag_scSPRITE2.txt");

// ---------------------------------------------------------------------------
// Realistic scSPRITE2 R2 sequence (125 bp)
//
// Layout: Y | SPACER | EVEN | SPACER | ODD | SPACER | EVEN | SPACER | ODD | SPACER | DPM
// Tags used (all exact matches, no laxity needed):
//   Y     NYbotLigEven_A1_Stg  TATTATGGT               (9 bp)
//   EVEN  Even2Bo1             ATACTGCGGCTGACG         (15 bp)
//   ODD   Odd2Bo1              TTCGTGGAATCTAGC         (15 bp)
//   EVEN  Even2Bo5             CTAGGTGGCGGTCTG         (15 bp)
//   ODD   Odd2Bo2              CCTACAGAAGTATCT         (15 bp)
//   DPM   DPM6bot1             TCATGTCTTCCGATCTTGGGTGTTTT (26 bp)
//   SPACER = 6 (default) → b"NNNNNN" between each tag pair
// ---------------------------------------------------------------------------
const REALISTIC_R2: &[u8] = b"\
TATTATGGT\
NNNNNN\
ATACTGCGGCTGACG\
NNNNNN\
TTCGTGGAATCTAGC\
NNNNNN\
CTAGGTGGCGGTCTG\
NNNNNN\
CCTACAGAAGTATCT\
NNNNNN\
TCATGTCTTCCGATCTTGGGTGTTTT";

// ---------------------------------------------------------------------------
// A: Hamming neighbor expansion
// ---------------------------------------------------------------------------

fn bench_generate_neighbors(c: &mut Criterion) {
    // Even2Bo1: representative 15-mer from the production config
    let seq = b"ATACTGCGGCTGACG";
    c.bench_function("generate_neighbors_15mer_k2", |b| {
        b.iter(|| generate_neighbors(std::hint::black_box(seq), 2))
    });
}

// ---------------------------------------------------------------------------
// B: check_read — full scSPRITE2 layout
// ---------------------------------------------------------------------------

fn bench_check_read(c: &mut Criterion) {
    let config   = Config::from_str_content(PROD_CONFIG).expect("PROD_CONFIG parse");
    let tag_maps = TagMaps::build(&config);

    c.bench_function("check_read_scsprite2_r2", |b| {
        b.iter(|| {
            let mut label = Vec::with_capacity(256);
            check_read(
                std::hint::black_box(REALISTIC_R2),
                &config.layout2,
                &tag_maps,
                config.spacer_len,
                config.laxity,
                &mut label,
            );
            label
        })
    });
}

// ---------------------------------------------------------------------------
// C: Batch throughput — 1 000 pairs processed in parallel via rayon
// ---------------------------------------------------------------------------

fn bench_batch_throughput(c: &mut Criterion) {
    const N: u64 = 1_000;

    let config   = Config::from_str_content(PROD_CONFIG).expect("PROD_CONFIG parse");
    let tag_maps = TagMaps::build(&config);

    let dummy_r1 = Record {
        header: b"@r1".to_vec(),
        seq:    b"ACGT".to_vec(),
        plus:   b"+".to_vec(),
        qual:   b"IIII".to_vec(),
    };
    let dummy_r2 = Record {
        header: b"@r1".to_vec(),
        seq:    REALISTIC_R2.to_vec(),
        plus:   b"+".to_vec(),
        qual:   vec![b'I'; REALISTIC_R2.len()],
    };
    let pairs: Vec<(Record, Record)> =
        (0..N).map(|_| (dummy_r1.clone(), dummy_r2.clone())).collect();

    let mut group = c.benchmark_group("batch_throughput");
    group.throughput(Throughput::Elements(N));
    group.bench_function(BenchmarkId::new("parallel", N), |b| {
        b.iter(|| {
            let _results: Vec<(Record, Record)> = std::hint::black_box(&pairs)
                .par_iter()
                .map(|(r1, r2)| process_pair(r1, r2, &config, &tag_maps))
                .collect();
        })
    });
    group.finish();
}

// ---------------------------------------------------------------------------
// D: Thread scaling — 10 000 pairs processed with 1 / 2 / 4 / 8 threads
//
// Uses non-global ThreadPool instances so each thread count is isolated from
// the others and from any global pool built by main().
// ---------------------------------------------------------------------------

fn bench_batch_scaling(c: &mut Criterion) {
    const N: u64 = 10_000;

    let config   = Config::from_str_content(PROD_CONFIG).expect("PROD_CONFIG parse");
    let tag_maps = TagMaps::build(&config);

    let dummy_r1 = Record {
        header: b"@r1".to_vec(),
        seq:    b"ACGT".to_vec(),
        plus:   b"+".to_vec(),
        qual:   b"IIII".to_vec(),
    };
    let dummy_r2 = Record {
        header: b"@r1".to_vec(),
        seq:    REALISTIC_R2.to_vec(),
        plus:   b"+".to_vec(),
        qual:   vec![b'I'; REALISTIC_R2.len()],
    };
    let pairs: Vec<(Record, Record)> =
        (0..N).map(|_| (dummy_r1.clone(), dummy_r2.clone())).collect();

    let mut group = c.benchmark_group("batch_scaling");
    group.throughput(Throughput::Elements(N));

    for &threads in &[1usize, 2, 4, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("ThreadPool build failed");

        group.bench_with_input(
            BenchmarkId::new("threads", threads),
            &threads,
            |b, _| {
                b.iter(|| {
                    pool.install(|| {
                        let _results: Vec<(Record, Record)> =
                            std::hint::black_box(&pairs)
                                .par_iter()
                                .map(|(r1, r2)| process_pair(r1, r2, &config, &tag_maps))
                                .collect();
                    })
                })
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------

criterion_group!(
    benches,
    bench_generate_neighbors,
    bench_check_read,
    bench_batch_throughput,
    bench_batch_scaling,
);
criterion_main!(benches);
