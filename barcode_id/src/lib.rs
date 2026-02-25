pub mod config;
pub mod fastq;
pub mod matcher;
pub mod tag;

#[cfg(test)]
mod parallel_tests {
    use std::collections::HashSet;
    use std::sync::{Arc, Mutex};

    use rayon::prelude::*;

    use crate::config::Config;
    use crate::fastq::Record;
    use crate::matcher::process_pair;
    use crate::tag::TagMaps;

    const PROD_CONFIG: &str =
        include_str!("../../misc/config_dpm6_y-stag_scSPRITE2.txt");

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

    /// Verify that process_pair is actually invoked from multiple distinct OS
    /// threads when rayon is given a 4-thread pool.  A pool of 4 threads
    /// processing 10 000 pairs must result in more than 1 unique thread ID
    /// being observed inside the closure — proving work was distributed.
    #[test]
    fn test_process_pair_runs_on_multiple_threads() {
        let config   = Config::from_str_content(PROD_CONFIG).expect("config parse");
        let tag_maps = TagMaps::build(&config);

        let r1 = Record { header: b"@r".to_vec(), seq: b"ACGT".to_vec(),
                          plus: b"+".to_vec(), qual: b"IIII".to_vec() };
        let r2 = Record { header: b"@r".to_vec(), seq: REALISTIC_R2.to_vec(),
                          plus: b"+".to_vec(), qual: vec![b'I'; REALISTIC_R2.len()] };
        let pairs: Vec<(Record, Record)> =
            (0..10_000).map(|_| (r1.clone(), r2.clone())).collect();

        let seen = Arc::new(Mutex::new(HashSet::new()));

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("ThreadPool build failed");

        pool.install(|| {
            pairs.par_iter().for_each(|(pr1, pr2)| {
                let _ = process_pair(pr1, pr2, &config, &tag_maps);
                seen.lock().unwrap().insert(std::thread::current().id());
            });
        });

        let n = seen.lock().unwrap().len();
        assert!(
            n > 1,
            "process_pair ran on only {} thread(s); expected >1 with a 4-thread pool",
            n,
        );
    }
}
