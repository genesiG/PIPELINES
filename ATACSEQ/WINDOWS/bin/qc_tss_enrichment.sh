#!/bin/bash
# qc_tss_enrichment.sh -- TSS enrichment of one normalised bigWig (deepTools).
# Called by TSS_ENRICHMENT (modules/local/tss_enrichment.nf):
#   qc_tss_enrichment.sh --sample S --bigwig S.<method>.<units>.bw --tss-bed tss.bed --threads N --flank 2000
# Arguments:
#   --sample      sample id
#   --bigwig      bigWig of the sample
#   --tss-bed     BED of transcription start sites (strand in column 6 is honoured by computeMatrix)
#   --threads     processors
#   --flank       bp on each side of the TSS (sets both --upstream and --downstream; default 2000).
#                 The ATAC pipeline passes tss_flank (legacy +-2 kb). NOT the same as --flank-bp.
#   --upstream    bp before the TSS (default: the --flank value; wins over --flank in either order)
#   --downstream  bp after the TSS (default: the --flank value; wins over --flank in either order)
#   --bin-size    matrix bin, bp (default 50; legacy binSize 50)
#   --center-bp   half-width of the "TSS" window, bp (default 250; bins whose centre lies within +-center-bp)
#   --flank-bp    width of each outer flank used as background, bp (default 200)
#   Other arguments go to computeMatrix (skipZeros and missingDataAsZero are always set, as in legacy).
# Steps: computeMatrix reference-point --referencePoint TSS --skipZeros --missingDataAsZero; plotProfile (SVG).
# Enrichment score, exact definition (mean profile = mean over regions of each bin, NaN values ignored; a bin with
#   no valid value counts as 0). Bin b has centre pos(b) = -upstream + (b + 0.5) * bin_bp, in bp from the TSS.
#     profile_max  = max over bins of the mean profile
#     profile_min  = min over bins of the mean profile
#     tss_signal   = max of the profile over bins with |pos| <= center-bp
#                    (profile_max if no bin falls in that window)
#     flank_mean   = mean of the profile over bins with pos < -upstream + flank-bp or pos > downstream - flank-bp
#     enrichment   = tss_signal / flank_mean                 (NA when flank_mean <= 0)
#     max_over_min = profile_max / profile_min               (NA when profile_min <= 0)
#   Legacy step_4qc_tss_enrichment.py score = max / min of the plotted mean profile, i.e. max_over_min.
#   enrichment is the template's centre-window / flank ratio; it has no legacy counterpart.
# Outputs: S.tss_enrichment.tsv (header + one row: sample_id n_regions window_up window_down bin_bp tss_signal
#   flank_mean enrichment profile_max profile_min max_over_min), S.tss_matrix.gz, S.tss_profile.svg
# Exit 3: empty matrix (no TSS covered by signal) or a bigWig/BED that is empty.
set -euo pipefail

sample="" bw="" tss="" threads=1 flank=2000 up="" down="" bin=50 center=250 fb=200 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)     sample="$2"; shift 2 ;;
        --bigwig)     bw="$2"; shift 2 ;;
        --tss-bed)    tss="$2"; shift 2 ;;
        --threads)    threads="$2"; shift 2 ;;
        --flank)      flank="$2"; shift 2 ;;
        --upstream)   up="$2"; shift 2 ;;
        --downstream) down="$2"; shift 2 ;;
        --bin-size)   bin="$2"; shift 2 ;;
        --center-bp)  center="$2"; shift 2 ;;
        --flank-bp)   fb="$2"; shift 2 ;;
        *)            extra+=("$1"); shift ;;
    esac
done
up="${up:-${flank}}"
down="${down:-${flank}}"
[[ -n "${sample}" && -s "${bw}" && -s "${tss}" ]] || { echo "qc_tss_enrichment.sh: missing --sample or empty --bigwig/--tss-bed" >&2; exit 3; }

computeMatrix reference-point \
    --referencePoint TSS \
    -R "${tss}" -S "${bw}" \
    -b "${up}" -a "${down}" --binSize "${bin}" \
    --skipZeros --missingDataAsZero \
    --numberOfProcessors "${threads}" \
    -o "${sample}.tss_matrix.gz" \
    ${extra[@]+"${extra[@]}"}

plotProfile -m "${sample}.tss_matrix.gz" \
    -out "${sample}.tss_profile.svg" --plotFileFormat svg \
    --samplesLabel "${sample}" --refPointLabel TSS \
    --plotTitle "${sample}: signal around TSS"

# matrix rows: chrom start end name score strand v1..vN ; line 1 is the "@{json}" header
gzip -dc "${sample}.tss_matrix.gz" | awk -v s="${sample}" -v up="${up}" -v down="${down}" -v bin="${bin}" -v cb="${center}" -v fb="${fb}" '
    BEGIN {OFS = "\t"}
    /^@/ {next}
    {
        n++
        for (i = 7; i <= NF; i++) {
            v = $i
            if (v == "nan" || v == "NaN" || v == "NA") continue
            sum[i - 7] += v; cnt[i - 7]++
        }
        nb = NF - 6
    }
    END {
        if (n == 0 || nb == 0) { print "empty TSS matrix: no TSS overlaps signal" > "/dev/stderr"; exit 3 }
        pmax = -1e300; pmin = 1e300; tmax = -1e300; fs = 0; fc = 0
        for (b = 0; b < nb; b++) {
            m = (cnt[b] > 0) ? sum[b] / cnt[b] : 0
            pos = -up + (b + 0.5) * bin
            if (m > pmax) pmax = m
            if (m < pmin) pmin = m
            ap = (pos < 0) ? -pos : pos
            if (ap <= cb && m > tmax) tmax = m
            if (pos < -up + fb || pos > down - fb) { fs += m; fc++ }
        }
        if (tmax < -1e299) tmax = pmax
        flankm = (fc > 0) ? fs / fc : 0
        enr = (flankm > 0) ? sprintf("%.4f", tmax / flankm) : "NA"
        mom = (pmin > 0) ? sprintf("%.4f", pmax / pmin) : "NA"
        print "sample_id", "n_regions", "window_up", "window_down", "bin_bp", "tss_signal", "flank_mean", "enrichment", "profile_max", "profile_min", "max_over_min"
        print s, n, up, down, bin, sprintf("%.6g", tmax), sprintf("%.6g", flankm), enr, sprintf("%.6g", pmax), sprintf("%.6g", pmin), mom
    }' > "${sample}.tss_enrichment.tsv"
