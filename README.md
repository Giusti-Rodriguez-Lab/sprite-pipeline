# scSPRITE Snakemake pipeline (UF Blue)

## How to run

```bash
sbatch sprite.sh
```

Monitor:
- Stdout: `snakemake.<JOBID>.out`
- Stderr: `snakemake.<JOBID>.err`

The batch script requests `--gres=gpu:1` for Parabricks and runs one Snakemake job at a time using up to 12 CPU cores (`-j 1 --cores 12`).

---

## Required inputs & folder layout

```
project/
├── Snakefile
├── sprite.sh
├── config.yaml
├── conda_envs.yaml
├── raw_fastq/
│   ├── SAMPLE_A_1.fq.gz
│   ├── SAMPLE_A_2.fq.gz
│   ├── SAMPLE_B_1.fq.gz
│   └── SAMPLE_B_2.fq.gz
```

**FASTQs**
- Place **paired-end** files in `raw_fastq/` named exactly:
  - `{sample}_1.fq.gz`
  - `{sample}_2.fq.gz`
- The pipeline auto-generates `raw_fastq/data_ID.txt` and builds the sample list from filenames.

**External tools / resources (via `config.yaml`)**
Provide full paths for:
- `RUNBC` — BarcodeID jar (`.jar`)
- `CONFIGBC` — BarcodeID config file (`.yaml`)
- `TRIMR1` — integer trim length for R1 (used by `awk substr`)
- `RUNLIG` — ligation efficiency script (`.py`)
- `RUNSAM` — `samtools` binary (or module path)
- `ALLTAG` — script to add barcode tags to BAM (`.py`)
- `MASK` — BED (repeats or regions to exclude)
- `RUNBED` — `bedtools` binary
- `GETCLUSTER` — clustering script (`.py`)
- `RUNQC` — `fastqc` binary
- `SPLITBYCELL` — script to split clusters per cell (`.py`)
- `NUMREADSCONTACTS` — script to compute per-cell reads/contacts (`.py`)
- `NORMCDF` — normalization script (`.py`)
- `CLUSTERBREAKDOWN` — script to summarize clusters (`.py`)
- `CHR` — list of chromosomes (used for printing)

**Genome / Parabricks**
- STAR genome-lib dir (default hardcoded):
  `/blue/.../Parabricks/STAR_index_100`
- FASTA:
  `/blue/.../GRCm38.primary_assembly.genome.fa`

---

## Minimal `config.yaml` template

```yaml
RUNBC: /path/to/BarcodeID.jar
CONFIGBC: /path/to/barcode_config.yaml
TRIMR1: 22
RUNLIG: /path/to/lig_efficiency.py
RUNSAM: /usr/bin/samtools
ALLTAG: /path/to/add_all_tags.py
MASK: /path/to/repeats.bed
RUNBED: /usr/bin/bedtools
GETCLUSTER: /path/to/make_clusters.py
RUNQC: /usr/bin/fastqc

SPLITBYCELL: /path/to/split_by_cell.py
NUMREADSCONTACTS: /path/to/num_reads_contacts.py
NORMCDF: /path/to/normcdf.py
CLUSTERBREAKDOWN: /path/to/cluster_breakdown.py

ASSEMBLY: mm10
RESOLUTION: 10000
DOWNWEIGHT: 1.0
RUNHIC: /path/to/hic_script.py
HEATMAP: /path/to/heatmap_script.py
MAX: 5
CELLSEL: 200
CHR: ["chr1","chr2","chrX","chrY"]
```

---

## What the pipeline produces (key outputs)

- **QC:**
  - `fastqc/` reports, and `log/fastqc_{sample}.txt`
- **Barcoding:**
  - `new_fastq/{sample}_R1.barcoded.fastq.gz`
  - `new_fastq/{sample}_R2.barcoded.fastq.gz`
- **Ligation efficiency:**
  - `ligation/{sample}.ligation_efficiency.txt`
  - `ligation/ligation_efficiency.txt`
- **Trimming (R1 only):**
  - `new_fastq/{sample}_R1.trimmed.fastq.gz`
- **Alignment (Parabricks → BAM):**
  - `alignment/{sample}Aligned.sortedByCoord.out.bam`
  - `alignment/{sample}Aligned.sortedByCoord.out.unique.bam`
  - `alignment/{sample}Aligned.sortedByCoord.out.unique.all_bcs.bam`
  - `alignment/{sample}Aligned.sortedByCoord.out.unique.all_bcs.masked.bam`
  - `alignment/all.bam`
- **Clustering:**
  - `cluster/clusters_all`
  - `cluster/clusters_all_reform`
  - `log/clusters_single_cells.txt`

---

## Notes / gotchas

- **Samples detection**: samples come from `raw_fastq/*_1.fq.gz`.
- **Conda**: rules that specify `conda: "conda_envs.yaml"` require that file.
- **Modules**: `sprite.sh` loads `java`, `snakemake/9.4.0`, `parabricks/4.5.1`, `conda/25.1.1`, and `samtools`.
- **Parabricks**: uses `pbrun rna_fq2bam` with GPU.
- **bedtools/samtools**: set paths in `config.yaml` or ensure binaries are on `$PATH`.
- **Hardcoded genome paths**: update `align_to_genome` rule if your reference differs.

Submit with:
```bash
sbatch sprite.sh
```
