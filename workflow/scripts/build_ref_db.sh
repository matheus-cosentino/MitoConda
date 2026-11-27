#!/usr/bin/env bash
############################################################################
#          Reference Mitochondrial Gene Database Extraction & Build        #
#                          MSc. Matheus Cosentino                          #
############################################################################
set -euo pipefail

ACTION="$1"
shift

if [[ "$ACTION" == "blast_hits" ]]; then
    GENE_QUERY="$1"
    TAXID="$2"
    OUT_TSV="$3"
    MAX_SP="$4"
    MAX_GEN="$5"
    MAX_GRP="$6"
    THREADS="$7"
    LOG="$8"

    export BLASTDB="resources/blast_db/mito:${BLASTDB:-}"
    > "${OUT_TSV}.raw"

    # Hierarchical lineage stratification across all major vertebrate classes, orders, suborders,
    # and families (all TaxIDs verified in resources/taxonomy/taxdump/names.dmp).
    # Each mitogenome hit climbs the NCBI taxonomy tree (Species -> Genus -> Subfamily -> Family ->
    # Suborder -> Order -> Class) and is assigned to its most specific matching lineage bucket so
    # hyper-abundant families never crowd out other orders or families.

    # 1. MAMMALIA (40674)
    # 1a. Chiroptera (9397) - all 20 bat families:
    CHIROPTERA_GROUPS="9398:20 58055:15 186994:15 9409:8 124754:6 217614:4 30562:12 59447:8 155036:4 94959:4 94957:6 59445:8 9415:30 148080:4 124757:5 59446:6 9436:18 981671:12 9431:28 3136025:4 9397:15"
    # 1b. Primates (9443) - Strepsirrhini, Tarsiiformes, Platyrrhini (Neotropical) & Catarrhini:
    PRIMATES_GROUPS="376915:18 376917:10 376912:6 9498:15 378855:12 376919:10 376918:8 9527:25 9577:10 9604:10 9443:10"
    # 1c. Rodentia (9989) - 33 suborders, families, and subfamilies:
    RODENTIA_GROUPS="10139:25 10161:35 10150:15 30646:15 723806:10 30645:25 10158:15 30649:15 108853:10 108856:10 10135:15 10167:20 74488:10 99599:10 10164:10 33550:20 40141:55 39087:35 337963:25 10026:20 337711:10 337677:15 39107:60 10045:25 326408:15 248797:10 10066:15 337664:25 337673:25 337669:10 687453:10 30648:30 1963761:15 55153:55 30650:20 51340:10 33553:10 29132:10 10015:25 10010:20 1963757:10 101663:10 10021:10 1963759:10 9989:25"
    # 1d. Lagomorpha (9975), Scandentia (9392) & Dermoptera (30656):
    EUARCHONTOGLIRES_OTHER="9979:20 9976:15 9975:10 9392:15 30656:6"
    # 1e. Eulipotyphla (9362) - Soricidae, Talpidae, Erinaceidae, Solenodontidae:
    EULIPOTYPHLA_GROUPS="9376:30 9373:18 9363:18 79803:6 9362:10"
    # 1f. Carnivora (33554) & Pholidota (9971):
    CARNIVORA_PHOLIDOTA_GROUPS="9681:22 9673:12 9697:12 379589:8 9676:6 379588:4 9608:20 9632:10 9655:24 9647:10 119825:8 379582:4 9702:10 9705:4 9709:12 33554:10 9972:12 9971:6"
    # 1g. Cetartiodactyla (91561) & Perissodactyla (9787):
    UNGULATA_CETACEA_GROUPS="9821:12 9827:8 9835:10 9831:5 9721:28 9846:8 30533:8 9850:22 9895:32 9892:5 9889:4 91561:10 9788:12 9799:8 9803:8 9787:6"
    # 1h. Afrotheria (311790), Xenarthra (9348), Marsupialia (9263) & Monotremata (9255):
    BASAL_MAMMALIA_GROUPS="9779:8 9774:6 9810:6 28734:15 9815:4 9369:15 9389:10 311790:10 9359:10 1811493:12 9352:8 227508:6 9349:10 9348:8 38605:25 38606:6 38607:4 38608:18 38611:10 38610:4 9307:18 9335:10 9331:8 37889:8 38624:4 9338:5 38609:10 9263:10 9256:4 9259:6 9255:4 40674:20"
    MAMMALIA_ALL="${CHIROPTERA_GROUPS} ${PRIMATES_GROUPS} ${RODENTIA_GROUPS} ${EUARCHONTOGLIRES_OTHER} ${EULIPOTYPHLA_GROUPS} ${CARNIVORA_PHOLIDOTA_GROUPS} ${UNGULATA_CETACEA_GROUPS} ${BASAL_MAMMALIA_GROUPS}"

    # 2. AVES (8782) - stratified across 24 major avian orders:
    AVES_GROUPS="8798:4 8792:4 8784:5 8819:5 8802:10 8976:25 8826:25 8929:20 8902:12 8892:18 8940:15 9108:18 8906:25 9230:12 30449:15 9205:18 3073854:12 2558200:25 30458:20 8936:15 9219:20 8948:15 9223:25 9126:90 8782:40"

    # 3. REPTILIA / NON-AVIAN SAUROPSIDA (Lepidosauria 8504, Testudines 8459, Crocodylia 1294634):
    LIZARDS_GROUPS="8560:22 8522:18 8530:15 88861:12 52185:12 81953:20 40248:15 8515:15 2024747:15 2024743:15 42424:12 8555:15 51964:12 8549:5 8504:25"
    SNAKES_GROUPS="34978:10 34977:8 8572:12 34984:12 8689:25 8602:25 8578:25 1390229:20 169862:15 1419188:12 8570:15"
    TESTUDINES_CROCODYLIA_GROUPS="8461:15 329146:8 44515:8 8487:18 328320:20 8476:20 8465:8 27792:3 34907:12 34912:10 8473:5 8459:10 1294636:10 8493:14 1294635:4 1294634:6"
    REPTILIA_ALL="${LIZARDS_GROUPS} ${SNAKES_GROUPS} ${TESTUDINES_CROCODYLIA_GROUPS}"

    # 4. AMPHIBIA (8292) - Anura (20 families), Caudata (8 families), Gymnophiona (5 families):
    ANURA_GROUPS="8352:10 30312:8 8343:8 8432:8 161696:18 8382:22 8418:25 43469:15 8372:20 613037:12 613053:8 507700:8 8427:18 8397:25 27577:20 43512:12 8412:10 685122:18 685128:8 685130:10 8342:15"
    CAUDATA_GYMNOPHIONA_GROUPS="30365:6 30366:18 8339:5 166053:6 8314:22 8294:12 8332:22 8310:4 8293:10 8446:10 8451:10 8454:6 1277737:8 30380:5 8445:8 8292:10"
    AMPHIBIA_ALL="${ANURA_GROUPS} ${CAUDATA_GYMNOPHIONA_GROUPS}"

    # 5. FISHES & BASAL VERTEBRATES (Actinopteri 186623, Chondrichthyes 7777, Sarcopterygii 8287, Cyclostomata 1476529):
    ACTINOPTERI_GROUPS="8288:6 7899:12 7921:4 41712:10 7933:15 32446:15 7952:25 7991:25 8002:12 7995:25 8006:15 8007:6 8043:12 1489878:15 1489911:20 28738:15 76071:10 8075:10 43697:8 1489900:10 1489907:12 8252:12 1489894:12 8111:20 31022:12 186623:25"
    CHONDRICHTHYES_BASAL_GROUPS="30483:15 30496:10 30503:8 282414:10 7858:10 117851:10 117888:6 117887:8 7864:6 7777:10 7894:3 2823314:6 7745:8 7761:6 7742:20"
    FISH_ALL="${ACTINOPTERI_GROUPS} ${CHONDRICHTHYES_BASAL_GROUPS}"

    if [[ "$TAXID" == "7742" ]]; then
        TARGET_GROUPS="${MAMMALIA_ALL} ${AVES_GROUPS} ${REPTILIA_ALL} ${AMPHIBIA_ALL} ${FISH_ALL}"
    elif [[ "$TAXID" == "32523" ]]; then
        TARGET_GROUPS="${MAMMALIA_ALL} ${AVES_GROUPS} ${REPTILIA_ALL} ${AMPHIBIA_ALL} 32523:25"
    elif [[ "$TAXID" == "40674" ]]; then
        TARGET_GROUPS="${MAMMALIA_ALL}"
    elif [[ "$TAXID" == "9989" ]]; then
        TARGET_GROUPS="${RODENTIA_GROUPS}"
    elif [[ "$TAXID" == "9397" ]]; then
        TARGET_GROUPS="${CHIROPTERA_GROUPS}"
    elif [[ "$TAXID" == "8782" ]]; then
        TARGET_GROUPS="${AVES_GROUPS}"
    elif [[ "$TAXID" == "8292" ]]; then
        TARGET_GROUPS="${AMPHIBIA_ALL}"
    elif [[ "$TAXID" == "8504" ]]; then
        TARGET_GROUPS="${LIZARDS_GROUPS} ${SNAKES_GROUPS}"
    else
        TARGET_GROUPS="${TAXID}:${MAX_GRP}"
    fi

    # Run a single high-capacity BLAST search against the target TaxID and include staxids
    blastn -task blastn \
           -word_size 11 \
           -query "$GENE_QUERY" \
           -db resources/blast_db/mito/mito \
           -taxids "$TAXID" \
           -outfmt "6 sseqid sstart send sstrand bitscore sscinames staxids" \
           -evalue 1e-5 \
           -max_target_seqs 100000 \
           -max_hsps 1 \
           -num_threads "$THREADS" > "${OUT_TSV}.raw" 2>> "$LOG" || true

    # Collect any reference sseqid already hit by sample BLAST searches so they are prioritized
    > "${OUT_TSV}.sample_hits"
    if compgen -G "output/*/individual_samples/*/Blast/*_Blastn_12s.txt" > /dev/null 2>&1; then
        awk -F'\t' '{print $2}' output/*/individual_samples/*/Blast/*_Blastn_12s.txt 2>/dev/null | sort -u > "${OUT_TSV}.sample_hits" || true
    elif compgen -G "output/*/samples/*/Blast/*_Blastn_12s.txt" > /dev/null 2>&1; then
        awk -F'\t' '{print $2}' output/*/samples/*/Blast/*_Blastn_12s.txt 2>/dev/null | sort -u > "${OUT_TSV}.sample_hits" || true
    fi

    # Climb the local NCBI taxonomy tree (taxonomy4blast.sqlite3 or nodes.dmp) bottom-up for each
    # BLAST hit's staxid to assign it to its most specific lineage bucket, and select balanced representatives
    python3 - "${OUT_TSV}.raw" "${OUT_TSV}.sample_hits" "$OUT_TSV" "$MAX_SP" "$MAX_GEN" "$TAXID" "$MAX_GRP" "$TARGET_GROUPS" << 'PYEOF'
import os
import sys
import sqlite3

raw_path, hits_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
max_sp, max_gen = int(sys.argv[4]), int(sys.argv[5])
root_taxid, default_cap = int(sys.argv[6]), int(sys.argv[7])
target_groups_str = sys.argv[8]

target_caps = {}
for entry in target_groups_str.split():
    if ":" in entry:
        tid_s, cap_s = entry.split(":", 1)
        target_caps[int(tid_s)] = int(cap_s)
if root_taxid not in target_caps:
    target_caps[root_taxid] = default_cap

sample_hits = set()
if os.path.exists(hits_path):
    with open(hits_path, "r") as f:
        for line in f:
            s = line.strip()
            if s:
                sample_hits.add(s)

# Parse raw BLAST hits and keep the best-scoring hit per sseqid
best_hits = {}
unique_staxids = set()
if os.path.exists(raw_path):
    with open(raw_path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            seqid, sstart, send, sstrand, bitscore_s, sciname, staxid_s = parts[:7]
            start, end = int(sstart), int(send)
            if start > end:
                start, end = end, start
            if (end - start) < 750:
                continue
            bitscore = float(bitscore_s)
            staxid = int(staxid_s.split(";")[0]) if staxid_s.split(";")[0].isdigit() else root_taxid
            if seqid not in best_hits or bitscore > best_hits[seqid][0]:
                best_hits[seqid] = (bitscore, seqid, start, end, sstrand, sciname, staxid)
                unique_staxids.add(staxid)

# Load taxonomy parent map from SQLite taxdb (or nodes.dmp fallback)
parent_map = {}
sqlite_path = "resources/blast_db/mito/taxonomy4blast.sqlite3"
nodes_path = "resources/taxonomy/taxdump/nodes.dmp"

if os.path.exists(sqlite_path):
    conn = sqlite3.connect(sqlite_path)
    cur = conn.cursor()
    pending = set(unique_staxids)
    visited = set()
    while pending:
        batch = list(pending)[:900]
        pending.difference_update(batch)
        visited.update(batch)
        placeholders = ",".join("?" for _ in batch)
        cur.execute(f"SELECT taxid, parent FROM TaxidInfo WHERE taxid IN ({placeholders})", batch)
        for tid, ptid in cur.fetchall():
            parent_map[int(tid)] = int(ptid)
            if int(ptid) not in visited and int(ptid) > 1:
                pending.add(int(ptid))
    conn.close()
elif os.path.exists(nodes_path):
    with open(nodes_path, "r") as f:
        for line in f:
            p = line.split("\t|\t", 2)
            if len(p) >= 2:
                parent_map[int(p[0])] = int(p[1])

bucket_cache = {}
def resolve_bucket(staxid):
    if staxid in bucket_cache:
        return bucket_cache[staxid]
    curr = staxid
    seen = set()
    while curr > 1 and curr not in seen:
        if curr in target_caps:
            bucket_cache[staxid] = curr
            return curr
        seen.add(curr)
        curr = parent_map.get(curr, 0)
    bucket_cache[staxid] = root_taxid
    return root_taxid

# Sort all unique mitogenomes by (is_sample_hit DESC, bitscore DESC)
sorted_records = sorted(
    best_hits.values(),
    key=lambda r: (1 if r[1] in sample_hits else 0, r[0]),
    reverse=True,
)

count_sp = {}
count_gen = {}
count_grp = {}

with open(out_path, "w") as out_f:
    for bitscore, seqid, start, end, sstrand, sciname, staxid in sorted_records:
        words = sciname.split()
        genus = words[0] if len(words) >= 1 else sciname
        species = f"{words[0]} {words[1]}" if len(words) >= 2 else sciname
        grp = resolve_bucket(staxid)
        grp_cap = target_caps.get(grp, default_cap)

        if seqid in sample_hits:
            if count_sp.get(species, 0) >= 1:
                continue
            count_sp[species] = count_sp.get(species, 0) + 1
            count_gen[genus] = count_gen.get(genus, 0) + 1
            count_grp[grp] = count_grp.get(grp, 0) + 1
            out_f.write(f"{seqid}\t{start}\t{end}\t{sstrand}\n")
            continue

        if max_sp > 0 and count_sp.get(species, 0) >= max_sp:
            continue
        if max_gen > 0 and count_gen.get(genus, 0) >= max_gen:
            continue
        if grp_cap > 0 and count_grp.get(grp, 0) >= grp_cap:
            continue

        count_sp[species] = count_sp.get(species, 0) + 1
        count_gen[genus] = count_gen.get(genus, 0) + 1
        count_grp[grp] = count_grp.get(grp, 0) + 1
        out_f.write(f"{seqid}\t{start}\t{end}\t{sstrand}\n")
PYEOF
    rm -f "${OUT_TSV}.raw" "${OUT_TSV}.sample_hits"

elif [[ "$ACTION" == "extract_raw" ]]; then
    BLAST_TSV="$1"
    RAW_FASTA="$2"
    LOG="$3"

    export BLASTDB="resources/blast_db/mito:${BLASTDB:-}"
    mkdir -p "$(dirname "$RAW_FASTA")"
    > "$RAW_FASTA"
    if [ -s "$BLAST_TSV" ]; then
        while IFS=$'\t' read -r seqid start end strand; do
            if [[ -z "$seqid" ]]; then continue; fi
            strand_flag="plus"
            [[ "$strand" == "minus" ]] && strand_flag="minus"
            if (( start > end )); then
                tmp="$start"
                start="$end"
                end="$tmp"
            fi
            blastdbcmd -db resources/blast_db/mito/mito \
                       -entry "$seqid" \
                       -range "$start-$end" \
                       -strand "$strand_flag" >> "$RAW_FASTA" 2>> "$LOG"
        done < "$BLAST_TSV"
    fi

elif [[ "$ACTION" == "normalize_headers" ]]; then
    RAW_FASTA="$1"
    GENE_FASTA="$2"

    mkdir -p "$(dirname "$GENE_FASTA")"
    awk 'BEGIN { OFS = "" }
         /^>/ {
             header = substr($0, 2)
             split(header, parts, " ")
             id = parts[1]
             gsub(/[:;=,\[\]()|-]/, "_", id)
             species = ""
             if (length(parts) >= 3) {
                 species = parts[2] "_" parts[3]
             } else if (length(parts) == 2) {
                 species = parts[2]
             }
             if (species == "") {
                 print ">" id
             } else {
                 print ">" id "_" species
             }
             next
         }
         { print }
    ' "$RAW_FASTA" > "$GENE_FASTA"
fi
