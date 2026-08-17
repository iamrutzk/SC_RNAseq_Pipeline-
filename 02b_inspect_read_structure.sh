#!/usr/bin/env bash
# Sanity check BEFORE running STARsolo: confirm which FASTQ is barcode/UMI
# (short read) vs cDNA (long read), and how long the barcode read actually
# is -- this determines --soloCBstart/len and --soloUMIstart/len.
# This matters a lot for inDrop data, where barcode structure is not a
# fixed 16+12 like 10x v3.
#
# Usage: 02b_inspect_read_structure.sh <R1.fastq.gz> <R2.fastq.gz>
set -euo pipefail

R1="$1"
R2="$2"

echo ">>> First 2 reads of $R1"
zcat "$R1" | head -n 8
echo
echo ">>> Read length distribution (first 1000 reads), $R1"
zcat "$R1" | awk 'NR%4==2{print length($0)}' | head -1000 | sort -n | uniq -c

echo
echo ">>> First 2 reads of $R2"
zcat "$R2" | head -n 8
echo
echo ">>> Read length distribution (first 1000 reads), $R2"
zcat "$R2" | awk 'NR%4==2{print length($0)}' | head -1000 | sort -n | uniq -c

cat <<'EOF'

>>> How to interpret:
  - Whichever file has short, FIXED-length reads (e.g. 16-28bp) is the
    barcode/UMI read -> feeds --soloCBwhitelist / --soloCBstart,len /
    --soloUMIstart,len.
  - Whichever file has long, cDNA-length reads (>50bp, variable-ish) is
    the transcript read -> the actual RNA-Seq input.
  - For 10x v2: barcode read = 26bp fixed (16 CB + 10 UMI).
  - For 10x v3: barcode read = 28bp fixed (16 CB + 12 UMI).
  - For inDrop: barcode read length and internal structure (CB1-linker-CB2-UMI)
    VARIES by library version. Do not assume 10x-style fixed CB/UMI split --
    consult the specific inDrop protocol used (check the GEO/SRA sample
    description and/or the original paper's methods) or use the "indrops"
    Python package (https://github.com/indrops/indrops) to pre-parse
    barcodes before feeding cDNA reads to STARsolo with
    --soloType CB_UMI_Simple and the correct, VERIFIED start/length values.
EOF
