# **Sequencing Metrics Summary — Documentation**

**Author:** Nathanael Herrera  
**Contact:** [ndh04c@gmail.com](mailto:ndh04c@gmail.com)

This document describes the metrics calculated by the **Qualimap + mosdepth + Picard summarizer** script (`summarize_seq_metrics_v5.py`). The script supports two modes: **Capture mode** for sequence capture/target enrichment data, and **WGS mode** for whole genome sequencing data.

---

## **Overview**

This pipeline pulls together information from three sources:

- **Qualimap** → read counts (mapped, unmapped, on-target, off-target)
- **mosdepth** → mean depth for target regions (Capture) or genome-wide (WGS)
- **Picard MarkDuplicates** → duplication rates and library complexity estimates

The script automatically determines which metrics to calculate based on whether you're running Capture or WGS mode.

---

## **Requirements**

### **Software Dependencies**

You'll need the following tools installed and available in your `$PATH`. The easiest way to set this up is through Conda:

```bash
# Use existing clean-map-qc conda environment (recommended)
# I just added mosdepth and matplotlib to my conda environment that has 
# includes raw read cleaning/ processing.

# Or you can create a new environment
conda create -n seq_metrics python=3.10
conda activate seq_metrics

# Install required tools
conda install -c bioconda qualimap mosdepth samtools bedtools
conda install conda-forge::matplotlib

# Assumes you have Picard tools
```

| Tool | Version Tested | Purpose |
|------|----------------|---------|
| **Python** | ≥3.10 | Script runtime |
| **mosdepth** | ≥0.3.3 | Coverage depth calculations |
| **samtools** | ≥1.15 | BAM file handling, read counting |
| **bedtools** | ≥2.30 | Fallback for near-target read counting |


### **Input File Requirements**

| Input | Description | How to Generate |
|-------|-------------|-----------------|
| **Qualimap results** | `genome_results.txt` files | Run `qualimap bamqc` on your BAMs |
| **BAM files** | Deduplicated, indexed BAMs | Your alignment pipeline |
| **Picard metrics** | `*_dedupe_metrics.txt` files | Run `picard MarkDuplicates` |
| **Reference FAI** | `.fai` index for your reference | `samtools faidx reference.fa` |
| **Target BED** | Capture targets (Capture mode only) | probe_coordinates_bed.bed |
| **Near-target BED** | Flanking regions (Capture mode, optional) | User-provided (see note below) |

**Note on near-target BED:** The `--near-bed` file should contain flanking regions around your targets (typically 100-250bp on each side) **with the target regions themselves excluded**. This is sometimes called a "Picard-style" near-target definition. How you generate this file is up to you. I used Bedtools to take the targets and add 200 bp to each end.

---

## **Usage**

### **Capture Mode**

For sequence capture / target enrichment data:

```bash
python3 summarize_seq_metrics_v5.py \
  -r /path/to/qualimap_results \
  -b targets.bed \
  --genome-fai reference.fa.fai \
  --bam-dir /path/to/bams \
  --picard-dir /path/to/picard_metrics \
  --cache-dir ./mosdepth_cache \
  --jobs 8 \
  -o capture_summary.csv
```

**Optional flags for Capture mode:**
- `--near-bed flanks.bed` — adds near-target metrics
- `--plot output.png` — generates a duplication analysis figure

### **WGS Mode**

For whole genome sequencing data:

```bash
python3 summarize_seq_metrics_v5.py \
  -r /path/to/qualimap_results \
  --wgs \
  --genome-fai reference.fa.fai \
  --bam-dir /path/to/bams \
  --picard-dir /path/to/picard_metrics \
  --cache-dir ./mosdepth_cache \
  --jobs 8 \
  -o wgs_summary.csv
```

You can also just omit the `--bed` flag entirely — the script will automatically switch to WGS mode.

### **All Command-Line Options**

| Flag | Required? | Description |
|------|-----------|-------------|
| `-r, --root` | Yes | Directory containing Qualimap `genome_results*.txt` files |
| `-b, --bed` | Capture only | Target BED file |
| `--wgs` | No | Force WGS mode (or just omit `--bed`) |
| `--genome-fai` | Yes | Reference `.fai` file |
| `--near-bed` | No | Flank-only BED for near-target metrics (Capture only) |
| `--bam-dir` | No | Directory with BAM files (if paths in Qualimap don't resolve) |
| `--picard-dir` | No | Directory with Picard MarkDuplicates output |
| `--cache-dir` | No | Where to store mosdepth results (default: `./mosdepth_summaries`) |
| `--jobs` | No | Parallel samples to process (default: 4) |
| `--timeout` | No | Seconds before near-target counting times out (default: 600) |
| `--force` | No | Re-run mosdepth even if cached results exist |
| `-o, --out` | No | Output CSV path (default: `sequencing_summary.csv`) |
| `--plot` | No | Generate duplication analysis plot (PNG path) |

---

## **Output Columns — Capture Mode**

| Column | What It Means | Units | Source | What to Look For |
|--------|---------------|-------|--------|------------------|
| **sample** | Sample identifier | String | From BAM filename | — |
| **Unmapped(reads)** | Reads that didn't align | Count | Qualimap | High = possible quality issues |
| **Total Mapped(reads)** | Reads aligned anywhere in genome | Count | Qualimap | Baseline for all percentages |
| **On-target(reads)** | Reads landing inside your targets | Count | Qualimap | Higher is better |
| **Off-target(reads)** | Reads landing outside your targets | Count | Mapped − On-target | Lower is better |
| **On-target rate (% of mapped reads)** | What fraction of mapped reads hit targets | Percent | (On-target ÷ Mapped) × 100 | **Key metric** — ≥70% is good |
| **Off-target rate (% of mapped reads)** | What fraction missed targets | Percent | (Off-target ÷ Mapped) × 100 | Complements on-target |
| **On-target (% of Total Reads)** | Overall efficiency including unmapped | Percent | (On-target ÷ Total) × 100 | Accounts for mapping rate |
| **Fold enrichment** | How much better than random sequencing | Fold (×) | On-target fraction ÷ (Target bp ÷ Genome bp) | ≥50× is good |
| **On-target(Mean Depth)** | Average coverage across targets | X | mosdepth | Higher = better variant calling |
| **Near-target(reads)** | Reads in flanking regions | Count | samtools/bedtools | Measures probe bleed-over |
| **Near-target rate (% of mapped reads)** | Fraction in flanks | Percent | (Near-target ÷ Mapped) × 100 | Some is expected |
| **Near-target(Mean Depth)** | Coverage in flanking regions | X | mosdepth | Usually lower than on-target |
| **Duplication rate** | PCR/optical duplicate fraction | Percent | Picard | <10% great, >30% not great |
| **Estimated Library Size** | Predicted unique molecules in library | Count | Picard | Higher = more complex library |
| **Optical Duplicates** | Duplicates from sequencer optics | Count | Picard | Should be low |
| **Unpaired Reads Examined** | Single-end reads checked for dups | Count | Picard | — |
| **Read Pairs Examined** | Paired-end reads checked for dups | Count | Picard | — |
| **Unique on-target reads** | Non-duplicate reads on targets | Count | On-target × (1 − dup rate) | **Usable data** |
| **qualimap_txt** | Path to source Qualimap file | Path | — | For troubleshooting |
| **mosdepth_summary** | Path to mosdepth output | Path | — | For troubleshooting |

---

## **Output Columns — WGS Mode**

WGS mode produces a simpler output since on-target/off-target concepts don't apply:

| Column | What It Means | Units | Source | What to Look For |
|--------|---------------|-------|--------|------------------|
| **sample** | Sample identifier | String | From BAM filename | — |
| **Unmapped(reads)** | Reads that didn't align | Count | Qualimap | High = possible quality issues |
| **Total Mapped(reads)** | Reads aligned to genome | Count | Qualimap | — |
| **Mapping rate (%)** | What fraction of reads aligned | Percent | (Mapped ÷ Total) × 100 | ≥95% typical for good data |
| **Mean Genome Depth** | Average coverage across genome | X | mosdepth | Higher is better |
| **Duplication rate** | PCR/optical duplicate fraction | Percent | Picard | <10% great, >30% not great |
| **Estimated Library Size** | Predicted unique molecules | Count | Picard | Higher = more complex library |
| **Optical Duplicates** | Duplicates from sequencer optics | Count | Picard | Should be low |
| **Unpaired Reads Examined** | Single-end reads checked | Count | Picard | — |
| **Read Pairs Examined** | Paired-end reads checked | Count | Picard | — |
| **Unique mapped reads** | Non-duplicate mapped reads | Count | Mapped × (1 − dup rate) | **Usable data** |
| **qualimap_txt** | Path to source file | Path | — | For troubleshooting |
| **mosdepth_summary** | Path to mosdepth output | Path | — | For troubleshooting |

---

## **Questions?**

Hit me up [ndh04c@gmail.com](mailto:ndh04c@gmail.com) if you run into issues or have suggestions for improvements.
