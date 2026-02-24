#!/usr/bin/env python3
"""Generate synthetic paired FASTQ.gz files for barcode_id benchmarking.

Reads the production barcode config to extract real tag sequences so that
all R2 reads contain exact barcode matches.  R1 is a fixed random sequence
(READ1 layout is absent in the production config, so its content is ignored).

Usage:
    python tests/gen_test_fastq.py \\
        --config  misc/config_dpm6_y-stag_scSPRITE2.txt \\
        --out-r1  /tmp/test_R1.fq.gz \\
        --out-r2  /tmp/test_R2.fq.gz \\
        --n-reads 500000
"""

import argparse
import gzip
import random
import sys


def parse_config(path):
    """Parse barcode config and return (layout2, spacer_len, first_tag_by_category).

    layout2     — list of category names from the READ2 = ... line
    spacer_len  — integer spacer length (default 6 if not in file)
    first_tag   — dict mapping CATEGORY -> first sequence seen for that category
    """
    layout2 = []
    spacer_len = 6      # DEFAULT_SPACER from config.rs
    first_tag = {}

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('READ2 = '):
                layout2 = [p.strip().upper() for p in line[len('READ2 = '):].split('|')]
            elif line.startswith('SPACER = '):
                spacer_len = int(line[len('SPACER = '):].strip())
            elif '\t' in line:
                fields = line.split('\t')
                if len(fields) >= 3:
                    cat = fields[0].strip().upper()
                    seq = fields[2].strip()
                    if cat not in first_tag:
                        first_tag[cat] = seq

    return layout2, spacer_len, first_tag


def build_r2_template(layout2, spacer_len, first_tag):
    """Build a single fixed R2 sequence from the layout (exact matches only)."""
    parts = []
    for cat in layout2:
        if cat == 'SPACER':
            parts.append('N' * spacer_len)
        else:
            seq = first_tag.get(cat)
            if seq is None:
                print(f"WARNING: no tag for category '{cat}', using 10× N", file=sys.stderr)
                seq = 'N' * 10
            parts.append(seq)
    return ''.join(parts)


def generate(out_r1, out_r2, r2_seq, n_reads, r1_len=50):
    r2_bytes = r2_seq.encode()
    r2_qual  = b'I' * len(r2_seq)
    r1_seq   = (''.join(random.choices('ACGT', k=r1_len))).encode()
    r1_qual  = b'I' * r1_len

    with gzip.open(out_r1, 'wb') as f1, gzip.open(out_r2, 'wb') as f2:
        for i in range(n_reads):
            name = f'read{i}'.encode()
            f1.write(b'@' + name + b'\n' + r1_seq + b'\n+\n' + r1_qual + b'\n')
            f2.write(b'@' + name + b'\n' + r2_bytes + b'\n+\n' + r2_qual + b'\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config',   required=True, help='Barcode config file')
    ap.add_argument('--out-r1',   required=True, help='Output R1 FASTQ.gz')
    ap.add_argument('--out-r2',   required=True, help='Output R2 FASTQ.gz')
    ap.add_argument('--n-reads',  type=int, default=500_000, help='Number of read pairs')
    ap.add_argument('--seed',     type=int, default=42,      help='Random seed')
    args = ap.parse_args()

    random.seed(args.seed)

    layout2, spacer_len, first_tag = parse_config(args.config)
    r2_seq = build_r2_template(layout2, spacer_len, first_tag)

    print(f"Layout:      {' | '.join(layout2)}", file=sys.stderr)
    print(f"Spacer len:  {spacer_len}", file=sys.stderr)
    print(f"R2 seq len:  {len(r2_seq)} bp", file=sys.stderr)
    print(f"R2 sequence: {r2_seq}", file=sys.stderr)
    print(f"Generating {args.n_reads:,} read pairs...", file=sys.stderr)

    generate(args.out_r1, args.out_r2, r2_seq, args.n_reads)

    print(f"Written: {args.out_r1}", file=sys.stderr)
    print(f"Written: {args.out_r2}", file=sys.stderr)


if __name__ == '__main__':
    main()
