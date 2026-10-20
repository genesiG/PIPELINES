#!/usr/bin/env python3
"""
cutadapt_stats.py
===================
Parses a cutadapt demultiplexing log (produced by DEMULTIPLEX_SE/PE) into
per-sample read-count and library-uniformity metrics (Gini, CV, min/max
fold). Ported from Scripts/step_1qc.py's run_cutadapt_stats() +
calc_uniformity_metrics(), made project-agnostic via a fixed samplesheet
schema (barcode_id,barcode,sample_name) instead of configurable column
indices.
"""

import argparse
import csv
import json
import math
import re


def calc_uniformity_metrics(counts: list[int]) -> dict:
    """Uniformity of per-sample read counts across demultiplexed samples.

    Zero-read samples (barcodes that got no reads) are included on purpose:
    an empty barcode is a demultiplexing failure and should lower uniformity.
    min_max_fold is None when the smallest sample has 0 reads (the fold is
    undefined; JSON null rather than inf, which json.dump would write as the
    non-standard token Infinity).
    """
    valid_counts = sorted([c for c in counts if c is not None and c >= 0])
    n = len(valid_counts)
    if n == 0 or sum(valid_counts) == 0:
        return {"gini": 0.0, "cv": 0.0, "min_max_fold": 0.0, "mean": 0.0, "median": 0.0}

    total = sum(valid_counts)
    mean = total / n
    stdev = math.sqrt(sum((c - mean) ** 2 for c in valid_counts) / n) if n > 1 else 0.0
    cv = stdev / mean if mean > 0 else 0.0

    gini = sum((2 * i - n - 1) * c for i, c in enumerate(valid_counts, 1)) / (n * total)

    min_val = valid_counts[0]
    max_val = valid_counts[-1]
    median = valid_counts[n // 2] if n % 2 != 0 else (valid_counts[n // 2 - 1] + valid_counts[n // 2]) / 2.0
    min_max_fold = (max_val / min_val) if min_val > 0 else None

    return {
        "gini": round(gini, 4),
        "cv": round(cv, 4),
        "min_max_fold": round(min_max_fold, 2) if min_max_fold is not None else None,
        "mean": round(mean, 1),
        "median": round(median, 1),
        "min": min_val,
        "max": max_val,
        "total": total,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log", required=True, help="cutadapt demux log")
    ap.add_argument("--samplesheet", required=True, help="samplesheet.csv (barcode_id,barcode,sample_name)")
    ap.add_argument("--project-name", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.log, "r", encoding="utf-8", errors="replace") as fh:
        log_text = fh.read()

    # ---- 1. Overall summary --------------------------------------------
    total_processed = 0
    reads_with_adapters = 0
    reads_discarded = 0
    total_bp_processed = 0
    total_bp_written = 0

    m = re.search(r"Total read pairs processed:\s+([0-9,]+)", log_text)
    if not m:
        m = re.search(r"Total reads processed:\s+([0-9,]+)", log_text)
    if m:
        total_processed = int(m.group(1).replace(",", ""))

    m = re.search(r"Read 1 with adapter:\s+([0-9,]+)", log_text)
    if not m:
        m = re.search(r"Reads with adapters:\s+([0-9,]+)", log_text)
    if m:
        reads_with_adapters = int(m.group(1).replace(",", ""))

    m = re.search(r"Pairs discarded as untrimmed:\s+([0-9,]+)", log_text)
    if not m:
        m = re.search(r"Reads discarded as untrimmed:\s+([0-9,]+)", log_text)
    if m:
        reads_discarded = int(m.group(1).replace(",", ""))

    m = re.search(r"Total basepairs processed:\s+([0-9,]+)\s+bp", log_text)
    if m:
        total_bp_processed = int(m.group(1).replace(",", ""))

    m = re.search(r"Total written \(filtered\):\s+([0-9,]+)\s+bp", log_text)
    if m:
        total_bp_written = int(m.group(1).replace(",", ""))

    pct_demux = round((reads_with_adapters / total_processed * 100.0), 2) if total_processed > 0 else 0.0
    pct_discarded = round((reads_discarded / total_processed * 100.0), 2) if total_processed > 0 else 0.0

    # ---- 2. Barcode/sample metadata & per-sample counts -------------------
    sample_barcodes = {}
    with open(args.samplesheet, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sname = (row.get("sample_name") or "").strip()
            if not sname:
                continue
            sample_barcodes[sname] = {
                "sequence": (row.get("barcode") or "").strip(),
                "id": (row.get("barcode_id") or sname).strip(),
            }

    sample_counts: dict[str, int] = {}

    # PE logs: "=== First read: Adapter NAME ===" ... "Trimmed: N times"
    pe_sections = re.findall(
        r"===\s+First read:\s+Adapter\s+([A-Za-z0-9_\-]+)\s+===.*?Trimmed:\s+([0-9,]+)\s+times",
        log_text, re.DOTALL,
    )
    for aname, cnt_str in pe_sections:
        sample_counts[aname] = int(cnt_str.replace(",", ""))

    if not pe_sections:
        se_sections = re.findall(
            r"===\s+Adapter\s+([A-Za-z0-9_\-]+)\s+===.*?Trimmed:\s+([0-9,]+)\s+times",
            log_text, re.DOTALL,
        )
        for aname, cnt_str in se_sections:
            sample_counts[aname] = int(cnt_str.replace(",", ""))

    sample_table = []
    total_trimmed_sum = sum(sample_counts.values())
    all_names = sorted(set(sample_barcodes.keys()) | set(sample_counts.keys()))

    for sname in all_names:
        cnt = sample_counts.get(sname, 0)
        pct_of_trimmed = round((cnt / total_trimmed_sum * 100.0), 2) if total_trimmed_sum > 0 else 0.0
        pct_of_total = round((cnt / total_processed * 100.0), 3) if total_processed > 0 else 0.0
        binfo = sample_barcodes.get(sname, {"sequence": "N/A", "id": "N/A"})
        sample_table.append({
            "sample_name": sname,
            "barcode_id": binfo["id"],
            "barcode_sequence": binfo["sequence"],
            "read_count": cnt,
            "pct_of_demuxed_reads": pct_of_trimmed,
            "pct_of_total_raw_reads": pct_of_total,
        })

    sample_table.sort(key=lambda x: x["read_count"], reverse=True)

    counts_list = [row["read_count"] for row in sample_table]
    uniformity = calc_uniformity_metrics(counts_list)

    metrics_data = {
        "project_name": args.project_name,
        "cutadapt_log_path": args.log,
        "overall_summary": {
            "total_reads_processed": total_processed,
            "reads_demultiplexed": reads_with_adapters,
            "pct_demultiplexed": pct_demux,
            "reads_discarded_untrimmed": reads_discarded,
            "pct_discarded": pct_discarded,
            "total_basepairs_processed": total_bp_processed,
            "total_basepairs_written": total_bp_written,
        },
        "sample_uniformity": uniformity,
        "samples": sample_table,
    }

    json_path = f"{args.outdir}/{args.project_name}_step1_demux_metrics.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(metrics_data, fh, indent=2)

    csv_path = f"{args.outdir}/{args.project_name}_step1_sample_counts.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "sample_name", "barcode_id", "barcode_sequence",
            "read_count", "pct_of_demuxed_reads", "pct_of_total_raw_reads",
        ])
        writer.writeheader()
        writer.writerows(sample_table)

    print(f"Saved: {json_path}")
    print(f"Saved: {csv_path}")
    print(
        f"Total Processed: {total_processed:,} | Demultiplexed: {reads_with_adapters:,} ({pct_demux}%) "
        f"| Discarded: {reads_discarded:,} ({pct_discarded}%)"
    )
    min_max_txt = "NA" if uniformity["min_max_fold"] is None else f"{uniformity['min_max_fold']}x"
    print(f"Uniformity Gini: {uniformity['gini']} | CV: {uniformity['cv']} | Min/Max fold: {min_max_txt}")


if __name__ == "__main__":
    main()
