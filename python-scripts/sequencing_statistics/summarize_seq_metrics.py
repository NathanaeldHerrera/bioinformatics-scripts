#!/usr/bin/env python3
"""
Summarize capture performance:
- Reads (Qualimap): total, mapped, on-target, off-target
- Rates as percentages: On-target, Off-target, Near-target (all % of mapped reads), plus On-target (% of total reads)
- Fold enrichment (on-target rate vs genome fraction targeted)
- On-target mean depth (mosdepth)
- Near-target mean depth (mosdepth on flank-only BED, Picard-style)
- Duplication rate (Picard MarkDuplicates) and Unique on-target reads

Author: Nathanael Herrera
Contact: ndh04c@gmail.com
"""

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Optional: plotting (only imported if --plot is used)
def import_plotting():
    """Lazily import matplotlib to avoid dependency if not plotting."""
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None


# ---------------------------- Qualimap parsing ---------------------------- #

def parse_genome_results(txt: Path):
    """
    Parse Qualimap genome_results.txt for read counts.
    Uses sections:
      >>>>>>> Globals
      >>>>>>> Globals inside
    """
    rec = {
        "sample": None,
        "bam": None,
        "reads_total": None,
        "reads_mapped": None,
        "reads_inside": None,
    }
    section = None
    for line in txt.read_text().splitlines():
        if line.startswith(">>>>>>>"):
            section = line.replace(">>>>>>>", "").strip().lower()
            continue

        m = re.search(r"bam file\s*=\s*(\S+)", line)
        if m and rec["bam"] is None:
            rec["bam"] = m.group(1)
            rec["sample"] = Path(rec["bam"]).stem.replace("_deduped", "").replace("_realigned", "")

        if section == "globals":
            m1 = re.search(r"number of reads\s*=\s*([\d,]+)", line, re.IGNORECASE)
            if m1 and rec["reads_total"] is None:
                rec["reads_total"] = int(m1.group(1).replace(",", ""))
            m2 = re.search(r"number of mapped reads\s*=\s*([\d,]+)", line, re.IGNORECASE)
            if m2 and rec["reads_mapped"] is None:
                rec["reads_mapped"] = int(m2.group(1).replace(",", ""))

        if section == "globals inside":
            m3 = re.search(r"number of mapped reads\s*=\s*([\d,]+)", line, re.IGNORECASE)
            if m3 and rec["reads_inside"] is None:
                rec["reads_inside"] = int(m3.group(1).replace(",", ""))

    return rec


# ---------------------------- mosdepth helpers ---------------------------- #

def parse_mosdepth_summary(summary_txt: Path):
    """
    Parse mosdepth summary supporting two formats:
      A) classic 2-line summary with labels 'total' and 'regions'
      B) per-chrom rows with header 'chrom length bases mean ...' and '*_region' rows for targets
    Returns: total_len, total_mean, reg_len, reg_mean (None if not derivable)
    """
    import re as _re
    if not summary_txt or not summary_txt.exists():
        return None, None, None, None

    with summary_txt.open() as fh:
        lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]

    if not lines:
        return None, None, None, None

    # Per-chrom format (with header)
    first = _re.split(r'[\t, ]+', lines[0])
    if first and first[0].lower() in {"chrom", "contig"}:
        rows = [_re.split(r'[\t, ]+', ln) for ln in lines[1:]]
        tot_len = tot_cov = 0
        reg_len = reg_cov = 0
        for parts in rows:
            if len(parts) < 4:
                continue
            label = parts[0]
            try:
                length = int(float(parts[1]))
                mean   = float(parts[3])
            except ValueError:
                continue
            if label.endswith("_region"):
                reg_len += length
                reg_cov += mean * length
            else:
                tot_len += length
                tot_cov += mean * length
        total_len  = tot_len if tot_len > 0 else None
        total_mean = (tot_cov / tot_len) if tot_len > 0 else None
        reg_len    = reg_len if reg_len > 0 else None
        reg_mean   = (reg_cov / reg_len) if reg_len and reg_len > 0 else None
        return total_len, total_mean, reg_len, reg_mean

    # Classic 'total'/'regions' labels (tabs/commas/spaces)
    total_len = total_mean = reg_len = reg_mean = None
    for ln in lines:
        parts = _re.split(r'[\t, ]+', ln)
        if len(parts) < 4:
            continue
        label = parts[0].lower()
        try:
            length = int(float(parts[1])); mean = float(parts[3])
        except ValueError:
            continue
        if label == "total":
            total_len, total_mean = length, mean
        elif label == "regions":
            reg_len, reg_mean = length, mean
    return total_len, total_mean, reg_len, reg_mean


def find_or_run_mosdepth(bam: str, bed: str, cache_dir: Path, sample_prefix: str, force: bool) -> Path | None:
    """
    Use cached <cache_dir>/<sample_prefix>.mosdepth.summary.txt if present (and not forced).
    Otherwise run mosdepth. If --fast-mode fails, retry without it.
    Deletes zero-byte summaries and prints stderr on failure.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = cache_dir / sample_prefix
    summary_path = Path(f"{out_prefix}.mosdepth.summary.txt")

    if summary_path.exists() and summary_path.stat().st_size == 0:
        summary_path.unlink(missing_ok=True)

    if summary_path.exists() and not force:
        return summary_path

    def run_cmd(cmd):
        print("[mosdepth]", " ".join(cmd))
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            if summary_path.exists() and summary_path.stat().st_size == 0:
                summary_path.unlink(missing_ok=True)
            print(f"[mosdepth][stderr for {sample_prefix}]\n{res.stderr.strip()}")
            raise subprocess.CalledProcessError(res.returncode, cmd, res.stdout, res.stderr)
        return True

    try:
        run_cmd(["mosdepth", "--by", bed, "--fast-mode", str(out_prefix), bam])
        if summary_path.exists() and summary_path.stat().st_size > 0:
            return summary_path
    except subprocess.CalledProcessError:
        print(f"[mosdepth] fast-mode failed for {sample_prefix}; retrying without --fast-mode...")

    try:
        run_cmd(["mosdepth", "--by", bed, str(out_prefix), bam])
        if summary_path.exists() and summary_path.stat().st_size > 0:
            return summary_path
    except subprocess.CalledProcessError:
        print(f"[mosdepth] FAILED for {sample_prefix}. Skipping.")
        return None

    if summary_path.exists() and summary_path.stat().st_size == 0:
        summary_path.unlink(missing_ok=True)
    return None


def fallback_regions_mean_from_bedgz(summary_path: Path):
    """
    Fallback: compute mean depth from <prefix>.regions.bed.gz
    Returns (reg_len, reg_mean) or (None, None) if not available.
    """
    import gzip
    if not summary_path:
        return None, None
    regions_bed_gz = Path(str(summary_path).replace(".mosdepth.summary.txt", ".regions.bed.gz"))
    if not regions_bed_gz.exists():
        return None, None
    tot_len = 0
    tot_cov = 0.0
    with gzip.open(regions_bed_gz, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                start = int(parts[1]); end = int(parts[2]); mean_depth = float(parts[3])
            except ValueError:
                continue
            length = max(0, end - start)
            tot_len += length
            tot_cov += mean_depth * length
    if tot_len > 0:
        return tot_len, (tot_cov / tot_len)
    return None, None


# ---------------------------- Picard metrics ----------------------------- #

def parse_picard_dup_rate(metrics_dir: Path, sample: str):
    """
    Parse Picard MarkDuplicates metrics to get duplication statistics.
    Expects a file named <sample>_dedupe_metrics.txt in metrics_dir.
    
    Returns a dict with:
        - dup_rate: PERCENT_DUPLICATION (fraction)
        - est_lib_size: ESTIMATED_LIBRARY_SIZE (may be None if not calculated)
        - unpaired_examined: UNPAIRED_READS_EXAMINED
        - pairs_examined: READ_PAIRS_EXAMINED
        - unpaired_dups: UNPAIRED_READ_DUPLICATES
        - pair_dups: READ_PAIR_DUPLICATES
        - optical_dups: READ_PAIR_OPTICAL_DUPLICATES
    Or (None, None) for backward compatibility if file not found.
    """
    cand = metrics_dir / f"{sample}_dedupe_metrics.txt"
    if not cand.exists():
        return None, None
    
    result = {
        'dup_rate': None,
        'est_lib_size': None,
        'unpaired_examined': None,
        'pairs_examined': None,
        'unpaired_dups': None,
        'pair_dups': None,
        'optical_dups': None,
    }
    
    with cand.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            if ln.startswith("LIBRARY"):
                header = ln.split("\t")
                for row in fh:
                    row = row.strip()
                    if not row or row.startswith("#"):
                        continue
                    vals = row.split("\t")
                    cols = {h: (vals[i] if i < len(vals) else "") for i, h in enumerate(header)}
                    
                    # Parse all available metrics
                    try:
                        result['dup_rate'] = float(cols.get("PERCENT_DUPLICATION", ""))
                    except ValueError:
                        pass
                    try:
                        est_str = cols.get("ESTIMATED_LIBRARY_SIZE", "").replace(",", "")
                        if est_str:
                            result['est_lib_size'] = int(est_str)
                    except ValueError:
                        pass
                    try:
                        result['unpaired_examined'] = int(cols.get("UNPAIRED_READS_EXAMINED", "").replace(",", ""))
                    except ValueError:
                        pass
                    try:
                        result['pairs_examined'] = int(cols.get("READ_PAIRS_EXAMINED", "").replace(",", ""))
                    except ValueError:
                        pass
                    try:
                        result['unpaired_dups'] = int(cols.get("UNPAIRED_READ_DUPLICATES", "").replace(",", ""))
                    except ValueError:
                        pass
                    try:
                        result['pair_dups'] = int(cols.get("READ_PAIR_DUPLICATES", "").replace(",", ""))
                    except ValueError:
                        pass
                    try:
                        result['optical_dups'] = int(cols.get("READ_PAIR_OPTICAL_DUPLICATES", "").replace(",", ""))
                    except ValueError:
                        pass
                    
                    # Return (dup_rate, est_lib_size) for backward compatibility,
                    # but also store full result in a way the caller can access
                    return result['dup_rate'], result['est_lib_size'], result
                break
    return None, None, None


def parse_picard_histogram(metrics_dir: Path, sample: str):
    """
    Parse the duplication histogram from Picard MarkDuplicates output.
    Returns list of (bin, all_sets, optical_sets, non_optical_sets) tuples.
    """
    cand = metrics_dir / f"{sample}_dedupe_metrics.txt"
    if not cand.exists():
        return None
    
    histogram = []
    in_histogram = False
    
    with cand.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if ln.startswith("## HISTOGRAM"):
                in_histogram = True
                continue
            if in_histogram:
                if not ln or ln.startswith("#"):
                    continue
                if ln.startswith("BIN"):
                    continue  # Skip header
                parts = ln.split("\t")
                if len(parts) >= 4:
                    try:
                        bin_val = float(parts[0])
                        all_sets = int(parts[2]) if len(parts) > 2 else 0
                        optical = int(parts[3]) if len(parts) > 3 else 0
                        non_optical = int(parts[4]) if len(parts) > 4 else 0
                        histogram.append((bin_val, all_sets, optical, non_optical))
                    except ValueError:
                        continue
    
    return histogram if histogram else None


# ---------------------------- BED / FAI utilities ------------------------ #

def load_chrom_sizes(fai_path: Path):
    sizes = {}
    with fai_path.open() as fh:
        for ln in fh:
            if not ln.strip():
                continue
            parts = ln.split()
            chrom, length = parts[0], int(parts[1])
            sizes[chrom] = length
    return sizes


def bed_union_length(bed_path: Path, chrom_sizes: dict | None = None):
    """
    Return total length of the UNION (merged, non-overlapping) of BED intervals.
    If chrom_sizes is given, clip intervals to contig bounds before merging.
    """
    by_chrom = {}
    with bed_path.open() as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            c = parts[0]
            try:
                s = int(parts[1]); e = int(parts[2])
            except ValueError:
                continue
            if chrom_sizes is not None:
                if c not in chrom_sizes:
                    continue
                L = chrom_sizes[c]
                if s < 0: s = 0
                if e > L: e = L
            if s >= e:
                continue
            by_chrom.setdefault(c, []).append((s, e))
    total = 0
    for c, ivals in by_chrom.items():
        ivals.sort()
        cur_s, cur_e = ivals[0]
        for s, e in ivals[1:]:
            if s <= cur_e:
                if e > cur_e:
                    cur_e = e
            else:
                total += (cur_e - cur_s)
                cur_s, cur_e = s, e
        total += (cur_e - cur_s)
    return total


# ---------------------------- Near-target (reads) ------------------------ #

def count_reads_in_bed(bam: Path, bed: Path, timeout: int = 600) -> int | None:
    """
    Count reads overlapping a BED file.
    
    Uses samtools view -c -L as primary method (faster, more memory efficient).
    Falls back to bedtools intersect -c if samtools fails.
    
    Args:
        bam: Path to BAM file
        bed: Path to BED file (e.g., near-target flanking regions)
        timeout: Maximum seconds to wait for command (default 600 = 10 min)
    
    Returns:
        Total read count overlapping BED regions, or None if all methods fail.
    
    Assumes BAM is deduplicated and indexed.
    """
    if not bed or not bed.exists():
        return None
    if not bam or not bam.exists():
        return None
    
    # Method 1: samtools view -c -L (fastest, most memory efficient)
    cmd_samtools = ["samtools", "view", "-c", "-L", str(bed), str(bam)]
    try:
        p = subprocess.run(cmd_samtools, capture_output=True, text=True, 
                          check=True, timeout=timeout)
        count = int(p.stdout.strip())
        return count
    except subprocess.TimeoutExpired:
        print(f"[near-target] samtools timed out after {timeout}s for {bam.name}, trying bedtools...")
    except subprocess.CalledProcessError as e:
        print(f"[near-target] samtools failed for {bam.name}: {e.stderr.strip()}, trying bedtools...")
    except ValueError as e:
        print(f"[near-target] samtools output parse error for {bam.name}: {e}, trying bedtools...")
    
    # Method 2: bedtools intersect -c (fallback, more memory but reliable)
    cmd_bedtools = ["bedtools", "intersect", "-c", "-a", str(bed), "-b", str(bam)]
    try:
        p = subprocess.run(cmd_bedtools, capture_output=True, text=True, 
                          check=True, timeout=timeout)
        total = 0
        for ln in p.stdout.splitlines():
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            try:
                cnt = int(parts[-1])
            except ValueError:
                continue
            total += cnt
        return total
    except subprocess.TimeoutExpired:
        print(f"[near-target] bedtools also timed out after {timeout}s for {bam.name}")
        return None
    except subprocess.CalledProcessError as e:
        print(f"[near-target] bedtools also failed for {bam.name}: {e.stderr.strip()}")
        return None


# ---------------------------- Per-sample processing ---------------------- #

def process_one(gr_path: Path,
                bed: Path,
                bam_dir: Path | None,
                cache_dir: Path,
                force: bool,
                picard_dir: Path | None,
                near_bed: Path | None,
                timeout: int = 600):
    """Process a single genome_results*.txt → summary row dict or None."""
    rec = parse_genome_results(gr_path)
    needed = [rec["bam"], rec["sample"], rec["reads_total"], rec["reads_mapped"], rec["reads_inside"]]
    if not all(needed):
        print(f"[warn] Skipping {gr_path} (missing required fields)")
        return None

    # Resolve BAM path via --bam-dir if provided
    bam_path = Path(rec["bam"])
    if bam_dir:
        bam_path = bam_dir / bam_path.name
    if not bam_path.exists():
        print(f"[warn] BAM not found for sample {rec['sample']}: {bam_path} — skipping")
        return None

    unmapped = rec["reads_total"] - rec["reads_mapped"]
    on_reads = rec["reads_inside"]
    off_reads = rec["reads_mapped"] - rec["reads_inside"]

    # Rates (fractions)
    on_rate = (on_reads / rec["reads_mapped"]) if rec["reads_mapped"] else None
    off_rate = (off_reads / rec["reads_mapped"]) if rec["reads_mapped"] else None
    pct_total_on = (on_reads / rec["reads_total"] * 100) if rec["reads_total"] else None

    # Duplication rate (Picard)
    dup_rate, est_lib_size, picard_metrics = (None, None, None)
    if picard_dir:
        dup_rate, est_lib_size, picard_metrics = parse_picard_dup_rate(picard_dir, rec["sample"])
    unique_on_reads = (on_reads * (1 - dup_rate)) if (dup_rate is not None) else None
    est_lib_size_val = est_lib_size  # Store for row output
    
    # Extract additional Picard metrics and calculate PCR vs optical duplication
    optical_dups = None
    unpaired_examined = None
    pairs_examined = None
    unpaired_dups = None
    pair_dups = None
    pcr_dup_rate = None
    optical_dup_rate = None
    
    if picard_metrics:
        optical_dups = picard_metrics.get('optical_dups')
        unpaired_examined = picard_metrics.get('unpaired_examined')
        pairs_examined = picard_metrics.get('pairs_examined')
        unpaired_dups = picard_metrics.get('unpaired_dups')
        pair_dups = picard_metrics.get('pair_dups')
        
        # Calculate optical and PCR duplicate rates
        # Total reads examined = unpaired + (pairs * 2)
        # Total duplicates = unpaired_dups + (pair_dups * 2)
        # Optical duplicates only counted for pairs (READ_PAIR_OPTICAL_DUPLICATES)
        if unpaired_examined is not None and pairs_examined is not None:
            total_examined = unpaired_examined + (pairs_examined * 2)
            
            if total_examined > 0:
                # Optical duplicate rate (optical dups * 2 because they're pairs)
                if optical_dups is not None:
                    optical_dup_rate = (optical_dups * 2) / total_examined
                
                # PCR duplicates = total duplicates - optical duplicates
                if unpaired_dups is not None and pair_dups is not None and optical_dups is not None:
                    total_dups = unpaired_dups + (pair_dups * 2)
                    optical_dups_as_reads = optical_dups * 2
                    pcr_dups = total_dups - optical_dups_as_reads
                    pcr_dup_rate = pcr_dups / total_examined

    # On-target mean depth (mosdepth on targets BED)
    on_md_summary = find_or_run_mosdepth(
        bam=str(bam_path), bed=str(bed), cache_dir=cache_dir, sample_prefix=rec["sample"], force=force
    )
    on_mean = None
    if on_md_summary and on_md_summary.exists():
        _, _, _, reg_mean = parse_mosdepth_summary(on_md_summary)
        if reg_mean is None:
            _, fallback_mean = fallback_regions_mean_from_bedgz(on_md_summary)
            on_mean = fallback_mean
        else:
            on_mean = reg_mean

    # Near-target: reads (via samtools/bedtools) & mean depth (mosdepth on flank-only BED)
    near_reads = None
    near_rate = None
    near_mean = None
    if near_bed and near_bed.exists():
        near_reads = count_reads_in_bed(bam_path, near_bed, timeout=timeout)
        if near_reads is not None and rec["reads_mapped"]:
            near_rate = near_reads / rec["reads_mapped"]

        near_md_summary = find_or_run_mosdepth(
            bam=str(bam_path), bed=str(near_bed), cache_dir=cache_dir, sample_prefix=f"{rec['sample']}_near", force=force
        )
        if near_md_summary and near_md_summary.exists():
            _, _, _, near_reg_mean = parse_mosdepth_summary(near_md_summary)
            if near_reg_mean is None:
                _, near_fallback_mean = fallback_regions_mean_from_bedgz(near_md_summary)
                near_mean = near_fallback_mean
            else:
                near_mean = near_reg_mean

    # Assemble row (format percentages here)
    row = {
        "sample": rec["sample"],
        "Unmapped(reads)": unmapped,
        "Total Mapped(reads)": rec["reads_mapped"],
        "On-target(reads)": on_reads,
        "Off-target(reads)": off_reads,
        "On-target rate (% of mapped reads)": f"{on_rate * 100:.2f}" if on_rate is not None else "",
        "Off-target rate (% of mapped reads)": f"{off_rate * 100:.2f}" if off_rate is not None else "",
        "On-target (% of Total Reads)": f"{pct_total_on:.2f}" if pct_total_on is not None else "",
        "Fold enrichment": "",  # filled later after genome/target lengths are known
        "On-target(Mean Depth)": f"{on_mean:.6f}" if on_mean is not None else "",
        "Near-target(reads)": f"{near_reads}" if near_reads is not None else "",
        "Near-target rate (% of mapped reads)": f"{near_rate * 100:.2f}" if near_rate is not None else "",
        "Near-target(Mean Depth)": f"{near_mean:.6f}" if near_mean is not None else "",
        "Duplication rate (total %)": f"{dup_rate * 100:.2f}" if dup_rate is not None else "",
        "PCR Dup rate (%)": f"{pcr_dup_rate * 100:.2f}" if pcr_dup_rate is not None else "",
        "Optical Dup rate (%)": f"{optical_dup_rate * 100:.2f}" if optical_dup_rate is not None else "",
        "Estimated Library Size": f"{est_lib_size_val}" if est_lib_size_val is not None else "",
        "Optical Duplicates": f"{optical_dups}" if optical_dups is not None else "",
        "Unpaired Reads Examined": f"{unpaired_examined}" if unpaired_examined is not None else "",
        "Read Pairs Examined": f"{pairs_examined}" if pairs_examined is not None else "",
        "Unique on-target reads": f"{unique_on_reads:.0f}" if unique_on_reads is not None else "",
        "qualimap_txt": str(gr_path),
        "mosdepth_summary": str(on_md_summary) if on_md_summary else "",
    }
    return row


# ----------------------------- Plotting -------------------------------- #

def generate_dup_vs_depth_plot(rows: list, output_path: str, plt):
    """
    Generate a scatter plot of duplication rate vs sequencing depth.
    
    Plots:
    1. Duplication rate (%) vs Total Mapped Reads (sequencing depth proxy)
    2. Duplication rate (%) vs Estimated Library Size (if available)
    3. Duplication rate (%) vs On-target Mean Depth
    
    Also fits a simple trend line to help visualize the relationship.
    """
    import numpy as np
    
    # Extract data
    samples = []
    dup_rates = []
    mapped_reads = []
    lib_sizes = []
    on_target_depths = []
    
    for r in rows:
        try:
            # Try both old and new column names for backward compatibility
            dup_str = r.get("Duplication rate (total %)", "") or r.get("Duplication rate", "") or ""
            dup = float(dup_str) if dup_str else float("nan")
            mapped = int(r.get("Total Mapped(reads)", 0))
            lib_size_str = r.get("Estimated Library Size", "")
            lib_size = int(lib_size_str) if lib_size_str else None
            depth_str = r.get("On-target(Mean Depth)", "")
            depth = float(depth_str) if depth_str else None
            
            if dup == dup and mapped > 0:  # NaN check
                samples.append(r.get("sample", ""))
                dup_rates.append(dup)
                mapped_reads.append(mapped)
                lib_sizes.append(lib_size)
                on_target_depths.append(depth)
        except (ValueError, TypeError):
            continue
    
    if len(dup_rates) < 2:
        print("[warn] Not enough data points for plotting")
        return
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Duplication Rate Analysis', fontsize=14, fontweight='bold')
    
    # Plot 1: Duplication vs Total Mapped Reads
    ax1 = axes[0, 0]
    ax1.scatter(np.array(mapped_reads) / 1e6, dup_rates, alpha=0.6, edgecolors='black', linewidth=0.5)
    ax1.set_xlabel('Total Mapped Reads (millions)')
    ax1.set_ylabel('Duplication Rate (%)')
    ax1.set_title('Duplication Rate vs Sequencing Depth')
    ax1.grid(True, alpha=0.3)
    
    # Add trend line
    x_arr = np.array(mapped_reads) / 1e6
    y_arr = np.array(dup_rates)
    if len(x_arr) > 1:
        z = np.polyfit(x_arr, y_arr, 1)
        p = np.poly1d(z)
        x_line = np.linspace(min(x_arr), max(x_arr), 100)
        ax1.plot(x_line, p(x_line), "r--", alpha=0.8, label=f'Trend: y={z[0]:.2f}x+{z[1]:.1f}')
        ax1.legend()
    
    # Calculate correlation
    corr = np.corrcoef(x_arr, y_arr)[0, 1]
    ax1.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax1.transAxes, 
             verticalalignment='top', fontsize=10, 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Plot 2: Duplication vs On-target Mean Depth
    ax2 = axes[0, 1]
    valid_depths = [(d, dup) for d, dup in zip(on_target_depths, dup_rates) if d is not None]
    if valid_depths:
        depths_arr = np.array([x[0] for x in valid_depths])
        dups_arr = np.array([x[1] for x in valid_depths])
        ax2.scatter(depths_arr, dups_arr, alpha=0.6, edgecolors='black', linewidth=0.5, color='green')
        ax2.set_xlabel('On-target Mean Depth (X)')
        ax2.set_ylabel('Duplication Rate (%)')
        ax2.set_title('Duplication Rate vs On-target Coverage')
        ax2.grid(True, alpha=0.3)
        
        if len(depths_arr) > 1:
            z2 = np.polyfit(depths_arr, dups_arr, 1)
            p2 = np.poly1d(z2)
            x_line2 = np.linspace(min(depths_arr), max(depths_arr), 100)
            ax2.plot(x_line2, p2(x_line2), "r--", alpha=0.8)
            corr2 = np.corrcoef(depths_arr, dups_arr)[0, 1]
            ax2.text(0.05, 0.95, f'r = {corr2:.3f}', transform=ax2.transAxes, 
                     verticalalignment='top', fontsize=10,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        ax2.text(0.5, 0.5, 'No depth data available', transform=ax2.transAxes,
                 ha='center', va='center')
    
    # Plot 3: Duplication vs Estimated Library Size
    ax3 = axes[1, 0]
    valid_libs = [(lib, dup, samp) for lib, dup, samp in zip(lib_sizes, dup_rates, samples) if lib is not None]
    if valid_libs:
        libs_arr = np.array([x[0] for x in valid_libs]) / 1e6
        dups_arr = np.array([x[1] for x in valid_libs])
        ax3.scatter(libs_arr, dups_arr, alpha=0.6, edgecolors='black', linewidth=0.5, color='purple')
        ax3.set_xlabel('Estimated Library Size (millions)')
        ax3.set_ylabel('Duplication Rate (%)')
        ax3.set_title('Duplication Rate vs Library Complexity')
        ax3.grid(True, alpha=0.3)
        
        if len(libs_arr) > 1:
            z3 = np.polyfit(libs_arr, dups_arr, 1)
            p3 = np.poly1d(z3)
            x_line3 = np.linspace(min(libs_arr), max(libs_arr), 100)
            ax3.plot(x_line3, p3(x_line3), "r--", alpha=0.8)
            corr3 = np.corrcoef(libs_arr, dups_arr)[0, 1]
            ax3.text(0.05, 0.95, f'r = {corr3:.3f}', transform=ax3.transAxes, 
                     verticalalignment='top', fontsize=10,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        ax3.text(0.5, 0.5, 'No library size data available', transform=ax3.transAxes,
                 ha='center', va='center')
    
    # Plot 4: Distribution of duplication rates (histogram)
    ax4 = axes[1, 1]
    ax4.hist(dup_rates, bins=20, edgecolor='black', alpha=0.7, color='steelblue')
    ax4.axvline(np.median(dup_rates), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(dup_rates):.1f}%')
    ax4.axvline(np.mean(dup_rates), color='orange', linestyle='--', linewidth=2, label=f'Mean: {np.mean(dup_rates):.1f}%')
    ax4.set_xlabel('Duplication Rate (%)')
    ax4.set_ylabel('Number of Samples')
    ax4.set_title('Distribution of Duplication Rates')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"[plot] Saved duplication analysis plot to {output_path}")
    plt.close()


# --------------------------------- Main ---------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description="Summarize capture enrichment metrics from Qualimap, mosdepth (targets & flank), and Picard."
    )
    ap.add_argument("-r", "--root", required=True,
                    help="Root dir to search recursively for genome_results*.txt (follows symlinks)")
    ap.add_argument("-b", "--bed", required=True,
                    help="Targets BED file used for mosdepth")
    ap.add_argument("--near-bed", default=None,
                    help="Flank-only BED (Picard-style near-target regions, excludes targets)")
    ap.add_argument("--genome-fai", required=True,
                    help="Reference .fai file (for genome size; also used to clip+union BED for fold enrichment)")
    ap.add_argument("--picard-dir", default=None,
                    help="Directory containing <sample>_dedupe_metrics.txt (Picard MarkDuplicates metrics)")
    ap.add_argument("-c", "--cache-dir", default="mosdepth_summaries",
                    help="Directory to write/reuse mosdepth summaries (default: mosdepth_summaries)")
    ap.add_argument("--bam-dir", default=None,
                    help="Directory containing final deduped BAMs (used to resolve BAM filenames)")
    ap.add_argument("--force", action="store_true",
                    help="Force re-running mosdepth even if summary exists")
    ap.add_argument("--jobs", type=int, default=4,
                    help="Number of samples to process in parallel (default: 4)")
    ap.add_argument("--timeout", type=int, default=600,
                    help="Timeout in seconds for near-target read counting per sample (default: 600)")
    ap.add_argument("-o", "--out", default="qualimap_targets_summary.csv",
                    help="Output CSV path")
    ap.add_argument("--plot", default=None,
                    help="Generate duplication vs sequencing depth plot (PNG output path)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    bed = Path(args.bed)
    near_bed = Path(args.near_bed).resolve() if args.near_bed else None
    cache_dir = Path(args.cache_dir)
    bam_dir = Path(args.bam_dir).resolve() if args.bam_dir else None
    picard_dir = Path(args.picard_dir).resolve() if args.picard_dir else None

    # Genome & target bp for enrichment
    fai = Path(args.genome_fai)
    sizes = load_chrom_sizes(fai)
    genome_bp = sum(sizes.values())
    target_bp = bed_union_length(bed, chrom_sizes=sizes)  # union length (clipped to FAI)

    # Gather Qualimap results (follow symlinks)
    qualimap_txts = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        for fname in filenames:
            if "genome_results" in fname and fname.endswith(".txt"):
                qualimap_txts.append(Path(dirpath) / fname)

    print(f"Found {len(qualimap_txts)} Qualimap result file(s) under {root}")
    for p in qualimap_txts[:10]:
        print("  -", p)
    if not qualimap_txts:
        print("No genome_results*.txt files found.")
        return

    # Parallel processing across samples
    rows = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futures = [
            ex.submit(process_one, gr, bed, bam_dir, cache_dir, args.force, picard_dir, near_bed, args.timeout)
            for gr in qualimap_txts
        ]
        for fut in as_completed(futures):
            row = fut.result()
            if row:
                rows.append(row)

    # Add enrichment columns (use FRACTION internally; the displayed on-target is percent)
    for r in rows:
        try:
            on_rate_pct = float(r.get("On-target rate (% of mapped reads)", "") or "nan")
            on_rate_frac = on_rate_pct / 100.0 if on_rate_pct == on_rate_pct else float("nan")  # NaN check
        except ValueError:
            on_rate_frac = float("nan")
        denom = (target_bp / genome_bp) if genome_bp > 0 else float("nan")
        fold = (on_rate_frac / denom) if (denom and denom > 0) else float("nan")
        r["Fold enrichment"] = f"{fold:.3f}" if (fold == fold) else ""  # NaN check

    # Write CSV
    cols = [
        "sample",
        "Unmapped(reads)", "Total Mapped(reads)", "On-target(reads)", "Off-target(reads)",
        "On-target rate (% of mapped reads)", "Off-target rate (% of mapped reads)",
        "On-target (% of Total Reads)", "Fold enrichment",
        "On-target(Mean Depth)",
        "Near-target(reads)", "Near-target rate (% of mapped reads)", "Near-target(Mean Depth)",
        "Duplication rate (total %)", "PCR Dup rate (%)", "Optical Dup rate (%)",
        "Estimated Library Size", "Optical Duplicates",
        "Unpaired Reads Examined", "Read Pairs Examined", "Unique on-target reads",
        "qualimap_txt", "mosdepth_summary",
    ]
    out_path = Path(args.out)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["sample"]):
            w.writerow({k: r.get(k, "") for k in cols})

    print(f"Wrote {len(rows)} rows to {out_path}")
    print(f"Genome bp (from FAI): {genome_bp:,}")
    print(f"Target union bp (from BED∩FAI): {target_bp:,}")

    # Generate duplication vs sequencing depth plot if requested
    if args.plot:
        plt = import_plotting()
        if plt is None:
            print("[warn] matplotlib not installed. Skipping plot generation.")
            print("       Install with: pip install matplotlib")
        else:
            generate_dup_vs_depth_plot(rows, args.plot, plt)


if __name__ == "__main__":
    main()

