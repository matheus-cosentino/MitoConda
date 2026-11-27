#!/usr/bin/env bash
############################################################################
#          Resolve Full Taxonomic Lineage for Reference Alignment          #
#                          MSc. Matheus Cosentino                          #
############################################################################
set -euo pipefail

REF_ALN="$1"
TAXDUMP_DIR="$2"
OUT_TSV="$3"

awk '
    /^>/ {
        hdr = substr($0, 2)
        split(hdr, a, " ")
        id = a[1]
        sp_part = id
        if (match(id, /^([A-Za-z0-9_]+[.][0-9]+)_[0-9]+_[0-9]+_(.+)$/, m)) {
            sp_part = m[2]
        }
        gsub(/^UNVERIFIED_/, "", sp_part)
        n = split(sp_part, w, "_")
        genus = w[1]
        species = (n >= 2) ? (w[1] " " w[2]) : w[1]
        print id "\t" species "\t" genus
    }
' "$REF_ALN" > "${OUT_TSV}.queries"

if [ -s "${OUT_TSV}.queries" ]; then
    taxonkit name2taxid -i 3 --data-dir "$TAXDUMP_DIR" "${OUT_TSV}.queries" | \
    awk -F'\t' 'BEGIN {OFS="\t"} $4 != "" && $4 + 0 > 0 {print $1, $2, $3, $4}' | \
    taxonkit lineage -i 4 --data-dir "$TAXDUMP_DIR" | \
    taxonkit reformat -i 5 --data-dir "$TAXDUMP_DIR" -f "{p};{c};{o};{f};{g}" > "${OUT_TSV}.resolved"

    awk -F'\t' '
        BEGIN { OFS="\t" }
        FNR == NR {
            full_lin = $5
            ref_lin = $6
            if (!($1 in best_ref) || (full_lin ~ /Chordata|Vertebrata/ && !(best_full[$1] ~ /Chordata|Vertebrata/))) {
                best_full[$1] = full_lin
                best_ref[$1] = ref_lin
            }
            next
        }
        {
            id = $1
            sp = $2
            gen = $3
            if (id in best_ref) {
                split(best_ref[id], t, ";")
                p = (t[1] != "") ? t[1] : "Unclassified"
                c = t[2]
                if (c == "" && best_full[id] ~ /Sauropsida/) c = "Sauropsida"
                if (c == "") c = "Unclassified"
                o = (t[3] != "") ? t[3] : "Unclassified"
                f = (t[4] != "") ? t[4] : "Unclassified"
                g = (t[5] != "") ? t[5] : gen
            } else {
                p = "Unclassified"
                c = "Unclassified"
                o = "Unclassified"
                f = "Unclassified"
                g = gen
            }
            print id, p, c, o, f, g, sp
        }
    ' "${OUT_TSV}.resolved" "${OUT_TSV}.queries" > "$OUT_TSV"
else
    > "$OUT_TSV"
fi
rm -f "${OUT_TSV}.queries" "${OUT_TSV}.resolved"

