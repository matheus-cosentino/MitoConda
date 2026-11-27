#!/usr/bin/env bash
############################################################################
#              Crop Reference-Only Columns, Run trimAl & Clean             #
#                          MSc. Matheus Cosentino                          #
############################################################################
set -euo pipefail

FRAGMENTS="$1"
ALIGNED="$2"
CLEANED="$3"
GT="$4"
LOG="$5"

if [ -s "$FRAGMENTS" ]; then
    awk '
        FNR == NR {
            if (/^>/) {
                hdr = substr($0, 2)
                split(hdr, a, " ")
                is_bc[a[1]] = 1
            }
            next
        }
        /^>/ {
            n++
            hdr = substr($0, 2)
            split(hdr, a, " ")
            ids[n] = a[1]
            headers[n] = $0
            seqs[n] = ""
            next
        }
        {
            seqs[n] = seqs[n] $0
        }
        END {
            if (n == 0) exit
            len = length(seqs[1])
            total_bc = 0
            for (i = 1; i <= n; i++) {
                if (is_bc[ids[i]]) {
                    total_bc++
                    for (p = 1; p <= len; p++) {
                        c = substr(seqs[i], p, 1)
                        if (c != "-" && c != ".") bc_cov[p]++
                    }
                }
            }
            if (total_bc == 0) {
                for (i = 1; i <= n; i++) {
                    print headers[i]
                    print seqs[i]
                }
            } else {
                k = 0
                for (p = 1; p <= len; p++) {
                    if (bc_cov[p] > 0) keep[++k] = p
                }
                for (i = 1; i <= n; i++) {
                    out_seq = ""
                    for (j = 1; j <= k; j++) {
                        out_seq = out_seq substr(seqs[i], keep[j], 1)
                    }
                    print headers[i]
                    print out_seq
                }
            }
        }
    ' "$FRAGMENTS" "$ALIGNED" > "${CLEANED}.cropped"
else
    cp "$ALIGNED" "${CLEANED}.cropped"
fi

trimal -in "${CLEANED}.cropped" -out "${CLEANED}.tmp" -gt "$GT" > "$LOG" 2>&1
if [ ! -s "${CLEANED}.tmp" ]; then
    cp "${CLEANED}.cropped" "${CLEANED}.tmp"
fi

# Enforce strict >=60% non-gap coverage across the trimmed alignment for all sequences
# so short fragments do not degrade phylogenetic alignment quality (dropped sample fragments
# are automatically inventoried and summarized in the MultiQC sequence filtering audit tables).
awk '
    /^>/ {
        if (n > 0) {
            seq_clean = seq
            gsub(/[-.Nn]/, "", seq_clean)
            if (length(seq) > 0 && (length(seq_clean) / length(seq)) >= 0.60) {
                print ">" cur_id
                print seq
            }
        }
        n++
        hdr = substr($0, 2)
        split(hdr, a, " ")
        id = a[1]
        if (match(id, /^([A-Za-z0-9_]+[.][0-9]+)_[0-9]+_[0-9]+_(.+)$/, m)) {
            id = m[2] "_" m[1]
        }
        if (seen[id]++) {
            id = id "_" seen[id]
        }
        cur_id = id
        seq = ""
        next
    }
    { seq = seq $0 }
    END {
        if (n > 0) {
            seq_clean = seq
            gsub(/[-.Nn]/, "", seq_clean)
            if (length(seq) > 0 && (length(seq_clean) / length(seq)) >= 0.60) {
                print ">" cur_id
                print seq
            }
        }
    }
' "${CLEANED}.tmp" > "$CLEANED"
rm -f "${CLEANED}.cropped" "${CLEANED}.tmp"

