############################################################################
#          Phylogenetic Patristic Distance & Richness Comparison           #
#                          MSc. Matheus Cosentino                          #
############################################################################

import os
import re
import math
import heapq
from collections import defaultdict
from Bio import Phylo, SeqIO


def _load_ref_taxonomy_lookup(ref_tax_path):
    """Load species and genus taxonomic lineage maps from ref_taxonomy.tsv."""
    genus_to_lineage = {}
    species_to_lineage = {}
    if ref_tax_path and os.path.exists(ref_tax_path):
        with open(ref_tax_path, "r") as rtf:
            for line in rtf:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 7:
                    _, p, c, o, f, g, s = parts[:7]
                    lin_tuple = (
                        p if p and p != "Unclassified" else "Unclassified",
                        c if c and c != "Unclassified" else "Unclassified",
                        o if o and o != "Unclassified" else "Unclassified",
                        f if f and f != "Unclassified" else "Unclassified",
                        g if g and g != "Unclassified" else "Unclassified",
                        s if s and s != "Unclassified" else "Unclassified",
                    )
                    if g and g != "Unclassified":
                        genus_to_lineage[g] = lin_tuple
                    if s and s != "Unclassified":
                        species_to_lineage[s] = lin_tuple
    return species_to_lineage, genus_to_lineage


def _infer_taxonomic_level_from_distance(p_dist, species, genus, family, order, cls):
    """
    Infer the most likely taxonomic relationship level between a sample cluster and its
    closest reference sequence based on empirical vertebrate 12S rRNA distance quartiles:
      - d <= 0.025 (<= 2.5% divergence, Intra-Species P75 = 0.022): Same Species
      - 0.025 < d <= 0.065 (Intra-Genus P75 = 0.062): Same Genus
      - 0.065 < d <= 0.125 (Intra-Family P75 = 0.123): Same Family
      - 0.125 < d <= 0.200 (Intra-Order range): Same Order (or Divergent / NUMT)
      - d > 0.200 (Intra-Class / Deep divergence): Same Class (Divergent / Putative NUMT)
    """
    g_label = genus if (genus and genus != "Unclassified") else (species.split()[0] if species else "Unknown")
    f_label = family if (family and family != "Unclassified") else g_label
    o_label = order if (order and order != "Unclassified") else f_label
    c_label = cls if (cls and cls != "Unclassified") else o_label

    if p_dist <= 0.025:
        return f"Same Species ({species})"
    elif p_dist <= 0.065:
        return f"Same Genus ({g_label} sp.)"
    elif p_dist <= 0.125:
        return f"Same Family ({f_label})"
    elif p_dist <= 0.200:
        return f"Same Order ({o_label} / Divergent)"
    else:
        return f"Same Class ({c_label} / Putative NUMT or Deeply Divergent)"


def run_phylo_distance_mqc(tree_path, ref_tax_path, dist_tsv, summary_mqc, species_mqc, bargraph_mqc, clusters_mqc, sample):
    """Calculate patristic genetic distance from each sample cluster to the nearest reference taxon."""
    os.makedirs(os.path.dirname(dist_tsv), exist_ok=True)

    species_to_lineage, genus_to_lineage = _load_ref_taxonomy_lookup(ref_tax_path)

    sample_prefix = f"{sample}_"
    short_sample = re.sub(r"_junto$", "", sample)
    short_prefix = f"{short_sample}_"
    tree = Phylo.read(tree_path, "newick")

    adj = defaultdict(list)
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            w = float(child.branch_length) if child.branch_length is not None else 0.0
            if w < 0.0:
                w = 0.0
            adj[id(parent)].append((id(child), w))
            adj[id(child)].append((id(parent), w))

    cluster_tips = []
    ref_tips = []
    for tip in tree.get_terminals():
        tip_name = tip.name or ""
        if tip_name.startswith(sample_prefix) or tip_name.startswith(short_prefix):
            cluster_tips.append(tip)
        elif not re.match(r"^barcode\d+_C\d+_size\d+", tip_name):
            ref_tips.append(tip)

    dist_to_ref = {}
    nearest_ref_name = {}
    pq = []
    for rtip in ref_tips:
        rid = id(rtip)
        rname = rtip.name or "Unknown_Ref"
        dist_to_ref[rid] = 0.0
        nearest_ref_name[rid] = rname
        heapq.heappush(pq, (0.0, rname, rid))

    while pq:
        d, rname, u = heapq.heappop(pq)
        if d > dist_to_ref.get(u, float("inf")):
            continue
        for v, w in adj[u]:
            nd = d + w
            if nd < dist_to_ref.get(v, float("inf")) or (
                nd == dist_to_ref.get(v, float("inf")) and rname < nearest_ref_name.get(v, "")
            ):
                dist_to_ref[v] = nd
                nearest_ref_name[v] = rname
                heapq.heappush(pq, (nd, rname, v))

    def parse_ref_label(ref_str):
        m_new = re.match(r"^(.+)_([A-Za-z]{2}_?\d+\.\d+)(?:_\d+)?$", ref_str)
        if m_new:
            species = m_new.group(1).replace("_", " ")
            acc = m_new.group(2)
            return acc, species
        m_old = re.match(r"^([A-Za-z0-9_]+\.\d+)_\d+_\d+_(.+)$", ref_str)
        if m_old:
            acc = m_old.group(1)
            species = m_old.group(2).replace("_", " ")
            return acc, species
        return ref_str, ref_str.replace("_", " ")

    def parse_cluster_info(tip_name):
        if tip_name.startswith(sample_prefix):
            raw_cluster = tip_name[len(sample_prefix):]
        elif tip_name.startswith(short_prefix):
            raw_cluster = tip_name[len(short_prefix):]
        else:
            raw_cluster = tip_name
        size_m = re.search(r"size[_=]?(\d+)", raw_cluster)
        size_val = int(size_m.group(1)) if size_m else 1
        c_m = re.search(r"^(C\d+_size\d+)", raw_cluster)
        uuid_m = re.search(r"centroid[_=]([0-9a-fA-F-]{8})", raw_cluster)
        if c_m:
            short_id = c_m.group(1)
        elif uuid_m:
            short_id = f"centroid_{uuid_m.group(1)}"
        else:
            short_id = raw_cluster[:24]
        return raw_cluster, short_id, size_val

    rows = []
    for ctip in cluster_tips:
        cid = id(ctip)
        tip_name = ctip.name or ""
        raw_cluster, short_id, size_val = parse_cluster_info(tip_name)
        p_dist = dist_to_ref.get(cid, 0.0)
        ref_str = nearest_ref_name.get(cid, "N/A")
        acc, species = parse_ref_label(ref_str)
        own_bl = float(ctip.branch_length) if ctip.branch_length is not None else 0.0

        gen_word = species.split()[0] if species else "Unclassified"
        lin = species_to_lineage.get(species) or genus_to_lineage.get(gen_word)
        if lin:
            _, c_val, o_val, f_val, g_val, _ = lin
        else:
            c_val, o_val, f_val, g_val = ("Unclassified", "Unclassified", "Unclassified", gen_word)

        lineage_str = f"{c_val} > {o_val} > {f_val} > {g_val}"
        inferred_level = _infer_taxonomic_level_from_distance(p_dist, species, g_val, f_val, o_val, c_val)

        rows.append({
            "sample": sample,
            "cluster_full": raw_cluster,
            "short_id": short_id,
            "reads": size_val,
            "closest_species": species,
            "closest_accession": acc,
            "closest_lineage": lineage_str,
            "putative_level": inferred_level,
            "c_val": c_val,
            "o_val": o_val,
            "f_val": f_val,
            "g_val": g_val,
            "closest_ref_raw": ref_str,
            "patristic_distance": p_dist,
            "terminal_branch_length": own_bl,
        })

    rows.sort(key=lambda r: (-r["reads"], r["patristic_distance"], r["short_id"]))

    with open(dist_tsv, "w") as f:
        f.write("Sample\tCluster_ID\tReads\tClosest_Species\tClosest_Accession\tClosest_Lineage\tPutative_Taxonomic_Level\tPatristic_Distance\tTerminal_Branch_Length\n")
        for r in rows:
            f.write(
                f"{r['sample']}\t{r['cluster_full']}\t{r['reads']}\t"
                f"{r['closest_species']}\t{r['closest_accession']}\t"
                f"{r['closest_lineage']}\t{r['putative_level']}\t"
                f"{r['patristic_distance']:.6f}\t{r['terminal_branch_length']:.6f}\n"
            )

    species_groups = defaultdict(list)
    for r in rows:
        species_groups[r["closest_species"]].append(r)

    n_clusters = len(rows)
    total_reads = sum(r["reads"] for r in rows)

    species_summaries = []
    for sp, sp_rows in species_groups.items():
        sp_reads = sum(x["reads"] for x in sp_rows)
        sp_clusters = len(sp_rows)
        sp_min = min(x["patristic_distance"] for x in sp_rows)
        sp_max = max(x["patristic_distance"] for x in sp_rows)
        sp_wmean = (
            sum(x["patristic_distance"] * x["reads"] for x in sp_rows) / sp_reads
            if sp_reads > 0 else 0.0
        )
        sp_pct = (sp_reads / total_reads * 100.0) if total_reads > 0 else 0.0
        accs = sorted(set(x["closest_accession"] for x in sp_rows))
        first_r = sp_rows[0]
        sp_putative = _infer_taxonomic_level_from_distance(
            sp_wmean, sp, first_r["g_val"], first_r["f_val"], first_r["o_val"], first_r["c_val"]
        )
        species_summaries.append({
            "species": sp,
            "accessions": ", ".join(accs),
            "lineage": first_r["closest_lineage"],
            "putative_level": sp_putative,
            "clusters": sp_clusters,
            "reads": sp_reads,
            "pct_reads": sp_pct,
            "min_dist": sp_min,
            "wmean_dist": sp_wmean,
            "max_dist": sp_max,
        })

    species_summaries.sort(key=lambda s: (-s["reads"], s["wmean_dist"], s["species"]))

    if n_clusters > 0:
        min_dist = min(r["patristic_distance"] for r in rows)
        max_dist = max(r["patristic_distance"] for r in rows)
        mean_dist = sum(r["patristic_distance"] for r in rows) / n_clusters
        weighted_dist = (
            sum(r["patristic_distance"] * r["reads"] for r in rows) / total_reads
            if total_reads > 0 else mean_dist
        )
        all_refs_str = " | ".join(
            f"{s['species']} ({s['reads']} reads, {s['pct_reads']:.1f}%, d={s['wmean_dist']:.4f})"
            for s in species_summaries
        )
    else:
        min_dist = max_dist = mean_dist = weighted_dist = 0.0
        all_refs_str = "None"

    with open(summary_mqc, "w") as f:
        f.write("# id: phylo_summary_stats\n")
        f.write("# section_name: '3.3 Phylogenetic Distance Summary (Per Sample)'\n")
        f.write(
            "# description: '<b>Methodology:</b> From the integrated FastTree Maximum-Likelihood phylogeny (GTR+GAMMA), a multi-source Dijkstra shortest-path algorithm computes the exact patristic genetic distance (sum of branch lengths in substitutions/site) from every VSEARCH consensus cluster to its closest NCBI reference leaf. "
            "<b>Columns &amp; Interpretation:</b> <b>Placed_Clusters</b> and <b>Cluster_Reads</b> show the total clusters and reads placed in the tree; "
            "<b>Identified_Ref_Taxa</b> is the number of distinct closest reference species; <b>All_Closest_References</b> lists every detected reference species with its read count, percentage of sample reads, and weighted mean patristic distance (<code>d</code>); "
            "<b>Min / Mean / Weighted_Mean / Max_Distance</b> summarize patristic distances across the sample (<code>d &le; 0.025</code> indicates a direct species-level match; <code>d = 0.025&ndash;0.125</code> indicates a close congeneric/confamilial relative not represented at the species level in NCBI; <code>d &gt; 0.125</code> typically indicates a nuclear mitochondrial pseudogene [NUMT] or a deeply divergent lineage).'\n"
        )
        f.write("# plot_type: 'table'\n")
        f.write("# pconfig:\n")
        f.write("#   id: 'phylo_summary_stats_table'\n")
        f.write("#   title: 'Phylogenetic Summary per Sample (FastTree)'\n")
        f.write("#   col1_header: 'Sample'\n")
        f.write("# headers:\n")
        f.write("#   Placed_Clusters:\n")
        f.write("#     title: 'Placed Clusters'\n")
        f.write("#     format: '{:,.0f}'\n")
        f.write("#   Cluster_Reads:\n")
        f.write("#     title: 'Cluster Reads'\n")
        f.write("#     format: '{:,.0f}'\n")
        f.write("#   Identified_Ref_Taxa:\n")
        f.write("#     title: 'Identified Ref Taxa'\n")
        f.write("#     format: '{:,.0f}'\n")
        f.write("#   All_Closest_References:\n")
        f.write("#     title: 'All Closest References'\n")
        f.write("#   Min_Distance:\n")
        f.write("#     title: 'Min Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("#   Mean_Distance:\n")
        f.write("#     title: 'Mean Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("#   Weighted_Mean_Distance:\n")
        f.write("#     title: 'Weighted Mean Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("#   Max_Distance:\n")
        f.write("#     title: 'Max Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("Sample\tPlaced_Clusters\tCluster_Reads\tIdentified_Ref_Taxa\tAll_Closest_References\tMin_Distance\tMean_Distance\tWeighted_Mean_Distance\tMax_Distance\n")
        f.write(
            f"{sample}\t{n_clusters}\t{total_reads}\t{len(species_summaries)}\t{all_refs_str}\t"
            f"{min_dist:.6f}\t{mean_dist:.6f}\t{weighted_dist:.6f}\t{max_dist:.6f}\n"
        )

    with open(species_mqc, "w") as f:
        f.write("# id: phylo_barcode_species_summary\n")
        f.write("# section_name: '3.4 Phylogenetic Placement: Barcodes vs. All Closest Reference Taxa'\n")
        f.write(
            "# description: '<b>Methodology:</b> Aggregates all consensus clusters within each sample that share the same closest phylogenetic reference species in the integrated FastTree ML tree. "
            "<b>Columns &amp; Interpretation:</b> Each row represents a unique <code>Sample | Closest_Species</code> pair, reporting its full taxonomic lineage (<code>Class &gt; Order &gt; Family &gt; Genus</code>), the inferred taxonomic confidence level (<b>Putative_Taxonomic_Level</b>) based on the read-weighted patristic distance, NCBI RefSeq accession(s), cluster count, read count, percentage of sample reads, and patristic genetic distances (substitutions/site).'\n"
        )
        f.write("# plot_type: 'table'\n")
        f.write("# pconfig:\n")
        f.write("#   id: 'phylo_barcode_species_summary_table'\n")
        f.write("#   title: 'Barcodes vs. Closest Reference Taxa'\n")
        f.write("#   col1_header: 'Barcode_Taxon'\n")
        f.write("# headers:\n")
        f.write("#   Sample:\n")
        f.write("#     title: 'Sample'\n")
        f.write("#   Closest_Species:\n")
        f.write("#     title: 'Closest Species'\n")
        f.write("#   Closest_Lineage:\n")
        f.write("#     title: 'Closest Lineage (Class > Order > Family > Genus)'\n")
        f.write("#   Putative_Taxonomic_Level:\n")
        f.write("#     title: 'Putative Taxonomic Level'\n")
        f.write("#   Reference_Accessions:\n")
        f.write("#     title: 'Reference Accessions'\n")
        f.write("#   Clusters:\n")
        f.write("#     title: 'Clusters'\n")
        f.write("#     format: '{:,.0f}'\n")
        f.write("#   Reads:\n")
        f.write("#     title: 'Reads'\n")
        f.write("#     format: '{:,.0f}'\n")
        f.write("#   Pct_Sample_Reads:\n")
        f.write("#     title: '% Sample Reads'\n")
        f.write("#     format: '{:,.2f}'\n")
        f.write("#   Min_Distance:\n")
        f.write("#     title: 'Min Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("#   Weighted_Mean_Distance:\n")
        f.write("#     title: 'Weighted Mean Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("#   Max_Distance:\n")
        f.write("#     title: 'Max Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("Barcode_Taxon\tSample\tClosest_Species\tClosest_Lineage\tPutative_Taxonomic_Level\tReference_Accessions\tClusters\tReads\tPct_Sample_Reads\tMin_Distance\tWeighted_Mean_Distance\tMax_Distance\n")
        for s in species_summaries:
            row_key = f"{sample} | {s['species']}"
            f.write(
                f"{row_key}\t{sample}\t{s['species']}\t{s['lineage']}\t{s['putative_level']}\t{s['accessions']}\t"
                f"{s['clusters']}\t{s['reads']}\t{s['pct_reads']:.2f}\t"
                f"{s['min_dist']:.6f}\t{s['wmean_dist']:.6f}\t{s['max_dist']:.6f}\n"
            )

    with open(bargraph_mqc, "w") as f:
        f.write("# id: phylo_closest_species_abundance\n")
        f.write("# section_name: '3.2 Phylogenetic Abundance by Closest Reference Species'\n")
        f.write(
            "# description: '<b>Methodology:</b> Stacked barplot of total read abundance per sample/barcode assigned to each closest reference species based on shortest patristic branch distance in the integrated FastTree Maximum-Likelihood tree. "
            "<b>How to Interpret:</b> Unlike BLAST+LCA (which collapses ambiguous top hits to higher ranks such as Family or Order), phylogenetic placement assigns every cluster to its single nearest evolutionary neighbor in the reference tree. "
            "Use the <i>Percentages</i> toggle in the plot toolbar to compare relative species composition across samples, and cross-check high-distance assignments in the tables below.'\n"
        )
        f.write("# plot_type: 'bargraph'\n")
        f.write("# pconfig:\n")
        f.write("#   id: 'phylo_closest_species_bargraph'\n")
        f.write("#   title: 'Reads per Closest Reference Species (Phylogeny)'\n")
        f.write("#   ylab: 'Number of Reads'\n")
        if species_summaries:
            sp_names = [s["species"] for s in species_summaries]
            sp_counts = [str(s["reads"]) for s in species_summaries]
            f.write("Sample\t" + "\t".join(sp_names) + "\n")
            f.write(f"{sample}\t" + "\t".join(sp_counts) + "\n")
        else:
            f.write("Sample\tUnplaced\n")
            f.write(f"{sample}\t0\n")

    with open(clusters_mqc, "w") as f:
        f.write("# id: phylo_cluster_distances\n")
        f.write("# section_name: '3.5 Phylogenetic Distance per Cluster (Closest Reference)'\n")
        f.write(
            "# description: '<b>Methodology:</b> Fine-grained cluster-level table reporting the phylogenetic placement, full taxonomic lineage, and distance-calibrated taxonomic assignment of every VSEARCH consensus cluster (<code>C01</code>, <code>C02</code>, ...) across all samples in the integrated tree. "
            "<b>Columns &amp; Taxonomic Thresholds (calibrated on the 12S Vertebrate Reference DB 75th-percentile bounds):</b> "
            "<b>Closest_Lineage</b> shows the NCBI hierarchy (<code>Class &gt; Order &gt; Family &gt; Genus</code>) of the nearest reference leaf; "
            "<b>Putative_Taxonomic_Level</b> classifies the likely taxonomic relationship from <b>Patristic_Distance</b> (<code>d</code>): "
            "<code>d &le; 0.025</code> &rarr; <b>Same Species</b>; "
            "<code>0.025 &lt; d &le; 0.065</code> &rarr; <b>Same Genus</b> (congeneric species not sequenced in DB); "
            "<code>0.065 &lt; d &le; 0.125</code> &rarr; <b>Same Family</b> (confamilial genus not in DB); "
            "<code>0.125 &lt; d &le; 0.200</code> &rarr; <b>Same Order</b> (divergent lineage within order); "
            "<code>d &gt; 0.200</code> &rarr; <b>Same Class (Putative NUMT or Deeply Divergent)</b>. "
            "<b>Terminal_Branch_Length</b> is the cluster&rsquo;s own private branch length.'\n"
        )
        f.write("# plot_type: 'table'\n")
        f.write("# pconfig:\n")
        f.write("#   id: 'phylo_cluster_distances_table'\n")
        f.write("#   title: 'Cluster Genetic Distance, Lineage & Putative Taxonomic Level'\n")
        f.write("#   col1_header: 'Cluster'\n")
        f.write("# headers:\n")
        f.write("#   Sample:\n")
        f.write("#     title: 'Sample'\n")
        f.write("#   Reads:\n")
        f.write("#     title: 'Reads'\n")
        f.write("#     format: '{:,.0f}'\n")
        f.write("#   Closest_Species:\n")
        f.write("#     title: 'Closest Species'\n")
        f.write("#   Closest_Accession:\n")
        f.write("#     title: 'Closest Accession'\n")
        f.write("#   Closest_Lineage:\n")
        f.write("#     title: 'Closest Lineage (Class > Order > Family > Genus)'\n")
        f.write("#   Putative_Taxonomic_Level:\n")
        f.write("#     title: 'Putative Taxonomic Level'\n")
        f.write("#   Patristic_Distance:\n")
        f.write("#     title: 'Patristic Distance'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("#   Terminal_Branch_Length:\n")
        f.write("#     title: 'Terminal Branch Length'\n")
        f.write("#     format: '{:,.4f}'\n")
        f.write("Cluster\tSample\tReads\tClosest_Species\tClosest_Accession\tClosest_Lineage\tPutative_Taxonomic_Level\tPatristic_Distance\tTerminal_Branch_Length\n")
        for r in rows:
            row_key = f"{r['sample']}_{r['short_id']}"
            f.write(
                f"{row_key}\t{r['sample']}\t{r['reads']}\t"
                f"{r['closest_species']}\t{r['closest_accession']}\t"
                f"{r['closest_lineage']}\t{r['putative_level']}\t"
                f"{r['patristic_distance']:.6f}\t{r['terminal_branch_length']:.6f}\n"
            )


def run_overall_abundance_richness_mqc(
    blast_files, lca_files, renamed_fastas, phylo_files, ref_tax_input,
    full_tsv, fasta_comp_mqc, richness_mqc, rank_abund_mqc, species_abund_mqc,
    samples, config, rank=None,
    filtered_audit_tsv=None, filter_summary_mqc=None, removed_seqs_mqc=None,
):
    """Generate per-FASTA method comparison, sequence filtering audit, and dataset-wide Abundance & Richness comparison tables (BLAST+LCA vs. Phylogeny)."""
    os.makedirs(os.path.dirname(full_tsv), exist_ok=True)

    rank_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
    rank_names = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
    if rank:
        defined_rank = str(rank).lower()
    else:
        cfg_ranks = config.get("taxonkit", {}).get("lca_ranks") or config.get("taxonkit", {}).get("lca_rank", ["order"])
        if isinstance(cfg_ranks, str):
            cfg_ranks = [cfg_ranks]
        defined_rank = str(cfg_ranks[0]).lower() if cfg_ranks else "order"
    if defined_rank not in rank_map:
        defined_rank = "order"
    rank_idx = rank_map[defined_rank]
    rank_label = defined_rank.title()

    ref_tax_path = ref_tax_input if isinstance(ref_tax_input, str) else (ref_tax_input[0] if ref_tax_input else "")
    species_to_lineage, genus_to_lineage = _load_ref_taxonomy_lookup(ref_tax_path)

    samples_list = list(samples)

    # 1. Parse BLAST top-hit info per (sample, qseqid)
    blast_top_hit = {}
    for s_name, b_path in zip(samples_list, blast_files):
        if not os.path.exists(b_path):
            continue
        with open(b_path, "r") as bf:
            for line in bf:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 14:
                    continue
                qid = parts[0]
                if (s_name, qid) in blast_top_hit:
                    continue
                sseqid = parts[1]
                acc_m = re.search(r"([A-Za-z]{2}_?\d+\.\d+)", sseqid)
                acc_str = acc_m.group(1) if acc_m else sseqid
                try:
                    pident_val = float(parts[2])
                except ValueError:
                    pident_val = 0.0
                sciname = parts[13].strip() if len(parts) > 13 and parts[13].strip() else acc_str
                blast_top_hit[(s_name, qid)] = (sciname, acc_str, pident_val)

    # 2. Parse renamed FASTA headers to map (sample, qseqid) <-> short cluster ID (C01_size2518)
    qid_to_short_cluster = {}
    short_cluster_to_qid = {}
    renamed_len_by_key = {}
    if renamed_fastas:
        for s_name, rf_path in zip(samples_list, renamed_fastas):
            if not os.path.exists(rf_path):
                continue
            short_s = re.sub(r"_junto$", "", s_name)
            for rec in SeqIO.parse(rf_path, "fasta"):
                # rec.id is e.g. barcode01_C01_size2518
                # rec.description is e.g. barcode01_junto|C01|size=2518|centroid=...;seqs=...;size=...
                c_short = rec.id
                if c_short.startswith(f"{short_s}_"):
                    c_short = c_short[len(short_s) + 1:]
                desc_parts = rec.description.split("|", 3)
                orig_qid = desc_parts[3].split()[0] if len(desc_parts) >= 4 else rec.id
                qid_to_short_cluster[(s_name, orig_qid)] = (rec.id, c_short)
                short_cluster_to_qid[(s_name, c_short)] = (rec.id, orig_qid)
                renamed_len_by_key[(s_name, orig_qid)] = len(rec.seq)

    # 3. Parse BLAST + LCA files per sample
    lca_sample_stats = {
        s: {
            "reads": 0,
            "clusters": 0,
            "rank_counts": defaultdict(int),
            "genera": set(),
            "species": set(),
        }
        for s in samples_list
    }
    lca_rank_agg = defaultdict(lambda: {"reads": 0, "clusters": 0, "samples": set(), "class": "Unclassified", "lineage": ""})
    lca_species_agg = defaultdict(lambda: {"reads": 0, "clusters": 0, "samples": set(), "lineage": ""})
    lca_records_by_key = {}

    for s_name, lca_path in zip(samples_list, lca_files):
        if not os.path.exists(lca_path):
            continue
        seen_q = set()
        with open(lca_path, "r") as lf:
            for line in lf:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                qid = parts[0]
                if qid in seen_q:
                    continue
                seen_q.add(qid)

                sz_m = re.search(r"size=(\d+)", qid)
                reads = int(sz_m.group(1)) if sz_m else 1

                full_lin = parts[-2] if len(parts) >= 2 else ""
                ref_lin = parts[-1]
                taxa = [t.strip() for t in ref_lin.split(";")]
                while len(taxa) < 6:
                    taxa.append("")
                if not taxa[1] and "Sauropsida" in full_lin:
                    taxa[1] = "Sauropsida"

                p_val = taxa[0] or "Unclassified"
                c_val = taxa[1] or "Unclassified"
                o_val = taxa[2] or "Unclassified"
                f_val = taxa[3] or "Unclassified"
                g_val = taxa[4] or ""
                s_val = taxa[5] or ""

                if g_val and g_val not in genus_to_lineage:
                    genus_to_lineage[g_val] = (p_val, c_val, o_val, f_val, g_val, s_val or g_val)

                non_empty = [t for t in taxa if t]
                if taxa[rank_idx]:
                    target_taxon = taxa[rank_idx]
                elif non_empty:
                    target_taxon = f"Unclassified ({non_empty[-1]})"
                else:
                    target_taxon = "Unclassified"

                # Determine lowest resolved rank in LCA
                lowest_lca_taxon = "Unclassified"
                lowest_lca_rank_name = "Unclassified"
                for r_i in range(5, -1, -1):
                    if taxa[r_i]:
                        lowest_lca_taxon = taxa[r_i]
                        lowest_lca_rank_name = rank_names[r_i]
                        break

                lowest_taxon = s_val if s_val else (f"{g_val} sp. (LCA)" if g_val else target_taxon)
                higher_str = f"{c_val} > {o_val} > {f_val}"
                full_lca_str = " > ".join(t for t in [c_val, o_val, f_val, g_val, s_val] if t and t != "Unclassified") or p_val

                lca_records_by_key[(s_name, qid)] = {
                    "sample": s_name,
                    "qid": qid,
                    "reads": reads,
                    "phylum": p_val,
                    "class": c_val,
                    "order": o_val,
                    "family": f_val,
                    "genus": g_val,
                    "species": s_val,
                    "defined_rank_taxon": target_taxon,
                    "lowest_lca_label": f"{lowest_lca_taxon} ({lowest_lca_rank_name})",
                    "lowest_lca_taxon": lowest_lca_taxon,
                    "lca_lineage": full_lca_str,
                }

                st = lca_sample_stats[s_name]
                st["reads"] += reads
                st["clusters"] += 1
                st["rank_counts"][target_taxon] += reads
                if g_val:
                    st["genera"].add(g_val)
                if s_val:
                    st["species"].add(s_val)

                ra = lca_rank_agg[target_taxon]
                ra["reads"] += reads
                ra["clusters"] += 1
                ra["samples"].add(s_name)
                if ra["class"] == "Unclassified" and c_val != "Unclassified":
                    ra["class"] = c_val
                if not ra["lineage"]:
                    ra["lineage"] = higher_str

                sa = lca_species_agg[lowest_taxon]
                sa["reads"] += reads
                sa["clusters"] += 1
                sa["samples"].add(s_name)
                if not sa["lineage"]:
                    sa["lineage"] = higher_str

    # 4. Parse Phylogeny genetic distance files per sample
    phylo_sample_stats = {
        s: {
            "reads": 0,
            "clusters": 0,
            "rank_counts": defaultdict(int),
            "species_counts": defaultdict(int),
            "genera": set(),
            "species": set(),
        }
        for s in samples_list
    }
    phylo_rank_agg = defaultdict(lambda: {
        "reads": 0, "clusters": 0, "samples": set(), "species_set": set(),
        "dist_weighted_sum": 0.0, "class": "Unclassified", "lineage": ""
    })
    phylo_species_agg = defaultdict(lambda: {
        "reads": 0, "clusters": 0, "samples": set(),
        "min_dist": float("inf"), "dist_weighted_sum": 0.0, "lineage": ""
    })
    phylo_records_by_short = {}

    if phylo_files:
        for s_name, p_path in zip(samples_list, phylo_files):
            if not os.path.exists(p_path):
                continue
            with open(p_path, "r") as pf:
                header_skipped = False
                for line in pf:
                    if not header_skipped:
                        header_skipped = True
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    try:
                        c_id = parts[1].strip()
                        reads = int(parts[2])
                        sp_name = parts[3].strip()
                        acc_name = parts[4].strip()
                        p_dist = float(parts[7]) if len(parts) >= 9 else float(parts[5])
                    except ValueError:
                        continue

                    gen_word = sp_name.split()[0] if sp_name else "Unclassified"
                    lin = species_to_lineage.get(sp_name) or genus_to_lineage.get(gen_word)
                    if lin:
                        p_val, c_val, o_val, f_val, g_val, _ = lin
                    else:
                        p_val, c_val, o_val, f_val, g_val = ("Chordata", "Unclassified", "Unclassified", "Unclassified", gen_word)

                    lin_by_idx = [p_val, c_val, o_val, f_val, g_val, sp_name]
                    target_taxon = lin_by_idx[rank_idx] if lin_by_idx[rank_idx] else "Unclassified"
                    higher_str = f"{c_val} > {o_val} > {f_val}"
                    full_phylo_lin = f"{c_val} > {o_val} > {f_val} > {g_val}"
                    putative_lvl = _infer_taxonomic_level_from_distance(p_dist, sp_name, g_val, f_val, o_val, c_val)

                    phylo_records_by_short[(s_name, c_id)] = {
                        "sample": s_name,
                        "c_id": c_id,
                        "reads": reads,
                        "species": sp_name,
                        "accession": acc_name,
                        "class": c_val,
                        "order": o_val,
                        "family": f_val,
                        "genus": g_val,
                        "lineage": full_phylo_lin,
                        "putative_level": putative_lvl,
                        "p_dist": p_dist,
                    }

                    pst = phylo_sample_stats[s_name]
                    pst["reads"] += reads
                    pst["clusters"] += 1
                    pst["rank_counts"][target_taxon] += reads
                    pst["species_counts"][sp_name] += reads
                    if g_val and g_val != "Unclassified":
                        pst["genera"].add(g_val)
                    if sp_name:
                        pst["species"].add(sp_name)

                    pra = phylo_rank_agg[target_taxon]
                    pra["reads"] += reads
                    pra["clusters"] += 1
                    pra["samples"].add(s_name)
                    pra["species_set"].add(sp_name)
                    pra["dist_weighted_sum"] += p_dist * reads
                    if pra["class"] == "Unclassified" and c_val != "Unclassified":
                        pra["class"] = c_val
                    if not pra["lineage"]:
                        pra["lineage"] = higher_str

                    psa = phylo_species_agg[sp_name]
                    psa["reads"] += reads
                    psa["clusters"] += 1
                    psa["samples"].add(s_name)
                    if p_dist < psa["min_dist"]:
                        psa["min_dist"] = p_dist
                    psa["dist_weighted_sum"] += p_dist * reads
                    if not psa["lineage"]:
                        psa["lineage"] = higher_str

    # 5. Write Section 4.1: Per-FASTA Sequence Comparison Table (BLAST Top Hit vs. TaxonKit LCA vs. Phylogeny)
    fasta_comp_rows = []
    for (s_name, qid), lca_rec in lca_records_by_key.items():
        short_s = re.sub(r"_junto$", "", s_name)
        b_sp, b_acc, b_pident = blast_top_hit.get((s_name, qid), ("N/A", "N/A", 0.0))
        b_top_label = f"{b_sp} ({b_acc})" if b_acc and b_acc != "N/A" else b_sp
        b_genus = b_sp.split()[0] if b_sp and b_sp != "N/A" else ""

        if (s_name, qid) in qid_to_short_cluster:
            fasta_id, c_short = qid_to_short_cluster[(s_name, qid)]
        else:
            uuid_m = re.search(r"centroid[_=]([0-9a-fA-F-]{8})", qid)
            short_u = uuid_m.group(1) if uuid_m else qid[:12]
            fasta_id = f"{short_s}_{short_u}_size{lca_rec['reads']}"
            c_short = None

        phylo_rec = phylo_records_by_short.get((s_name, c_short)) if c_short else None
        if phylo_rec:
            p_sp = phylo_rec["species"]
            p_acc = phylo_rec["accession"]
            p_sp_label = f"{p_sp} ({p_acc})"
            p_dist = phylo_rec["p_dist"]
            p_putative = phylo_rec["putative_level"]
            p_genus = phylo_rec["genus"]
            p_family = phylo_rec["family"]
            p_order = phylo_rec["order"]

            # Evaluate concordance across BLAST Top Hit, LCA, and Phylogeny
            if p_dist > 0.125:
                concordance = f"Divergent / Putative NUMT (d={p_dist:.4f}, nearest: {p_sp})"
            elif b_sp == p_sp and p_dist <= 0.025:
                concordance = f"Exact Species Match ({p_sp})"
            elif (b_genus and b_genus == p_genus) or (lca_rec["genus"] and lca_rec["genus"] == p_genus):
                concordance = f"Genus Concordant ({p_genus}; BLAST={b_sp}, Phylo={p_sp})"
            elif lca_rec["family"] and lca_rec["family"] == p_family:
                concordance = f"Family Concordant ({p_family}; LCA={lca_rec['lowest_lca_taxon']}, Phylo={p_sp})"
            elif lca_rec["order"] and lca_rec["order"] == p_order:
                concordance = f"Order Concordant ({p_order}; LCA={lca_rec['lowest_lca_taxon']}, Phylo={p_sp})"
            else:
                concordance = f"Discrepant (BLAST={b_sp}, Phylo={p_sp})"
        else:
            p_sp_label = "Not Placed in Tree (<60% Alignment Overlap or Non-Target Clade)"
            p_dist = 0.0
            p_putative = "BLAST + LCA Only"
            if c_short is not None:
                concordance = "BLAST + LCA Only (Filtered by <60% Alignment Coverage)"
            else:
                concordance = f"BLAST + LCA Only (Outside Target Clade: {lca_rec['class']})"

        fasta_comp_rows.append({
            "fasta_id": fasta_id,
            "sample": s_name,
            "reads": lca_rec["reads"],
            "blast_top_hit": b_top_label,
            "blast_pident": b_pident,
            "lca_assigned": lca_rec["lowest_lca_label"],
            "lca_lineage": lca_rec["lca_lineage"],
            "phylo_closest": p_sp_label,
            "phylo_putative": p_putative,
            "phylo_dist": p_dist,
            "concordance": concordance,
        })

    fasta_comp_rows.sort(key=lambda r: (r["sample"], -r["reads"], r["fasta_id"]))

    with open(fasta_comp_mqc, "w") as fcf:
        fcf.write("# id: fasta_method_comparison\n")
        fcf.write("# section_name: '4.1 Per-FASTA Sequence Comparison (BLAST Top Hit vs. LCA vs. Phylogeny)'\n")
        fcf.write(
            "# description: '<b>Methodology:</b> Direct sequence-by-sequence comparison of every VSEARCH 12S consensus FASTA cluster across all three identification layers: "
            "<b>(1) BLAST Top Hit</b> (highest-scoring mitochondrial match and its % sequence identity), "
            "<b>(2) TaxonKit LCA</b> (conservative lowest common ancestor across all top BLAST hits), and "
            "<b>(3) FastTree Phylogenetic Placement</b> (closest reference species in the integrated ML tree, patristic distance <code>d</code>, and distance-calibrated taxonomic assignment). "
            "<b>How to Interpret:</b> Look at <b>Method_Concordance</b> to verify whether a FASTA sequence is an <i>Exact Species Match</i> across methods, a <i>Genus/Family Concordant</i> sequence where Phylogeny resolves a conservative family-level LCA call, a <i>Divergent / Putative NUMT</i> co-amplification (<code>d &gt; 0.125</code>), or a short fragment retained only in <i>BLAST+LCA</i>.'\n"
        )
        fcf.write("# plot_type: 'table'\n")
        fcf.write("# pconfig:\n")
        fcf.write("#   id: 'fasta_method_comparison_table'\n")
        fcf.write("#   title: 'Per-FASTA 12S Sequence Comparison: BLAST Top Hit vs. LCA vs. Phylogeny'\n")
        fcf.write("#   col1_header: 'FASTA_Sequence_ID'\n")
        fcf.write("# headers:\n")
        fcf.write("#   Sample:\n")
        fcf.write("#     title: 'Sample'\n")
        fcf.write("#   Reads:\n")
        fcf.write("#     title: 'Reads'\n")
        fcf.write("#     format: '{:,.0f}'\n")
        fcf.write("#   BLAST_Top_Hit:\n")
        fcf.write("#     title: 'BLAST Top Hit (Accession)'\n")
        fcf.write("#   BLAST_Identity_Pct:\n")
        fcf.write("#     title: 'BLAST Identity (%)'\n")
        fcf.write("#     format: '{:,.2f}'\n")
        fcf.write("#   LCA_Assigned_Taxon:\n")
        fcf.write("#     title: 'LCA Consensus Taxon (Rank)'\n")
        fcf.write("#   LCA_Lineage:\n")
        fcf.write("#     title: 'LCA Lineage'\n")
        fcf.write("#   Phylo_Closest_Species:\n")
        fcf.write("#     title: 'Phylogeny Closest Species'\n")
        fcf.write("#   Phylo_Putative_Level:\n")
        fcf.write("#     title: 'Phylogeny Putative Level'\n")
        fcf.write("#   Phylo_Patristic_Distance:\n")
        fcf.write("#     title: 'Patristic Distance'\n")
        fcf.write("#     format: '{:,.4f}'\n")
        fcf.write("#   Method_Concordance:\n")
        fcf.write("#     title: 'Method Concordance (BLAST + LCA vs. Phylogeny)'\n")
        fcf.write(
            "FASTA_Sequence_ID\tSample\tReads\tBLAST_Top_Hit\tBLAST_Identity_Pct\tLCA_Assigned_Taxon\tLCA_Lineage\t"
            "Phylo_Closest_Species\tPhylo_Putative_Level\tPhylo_Patristic_Distance\tMethod_Concordance\n"
        )
        for r in fasta_comp_rows:
            fcf.write(
                f"{r['fasta_id']}\t{r['sample']}\t{r['reads']}\t{r['blast_top_hit']}\t{r['blast_pident']:.2f}\t"
                f"{r['lca_assigned']}\t{r['lca_lineage']}\t{r['phylo_closest']}\t{r['phylo_putative']}\t"
                f"{r['phylo_dist']:.6f}\t{r['concordance']}\n"
            )

    # 5b. Sequence Filtering & Retention Audit (VSEARCH Consensus -> BLAST -> Target Clade -> trimAl / Phylogeny)
    if filtered_audit_tsv or filter_summary_mqc or removed_seqs_mqc:
        audit_summary_per_sample = []
        all_audit_seq_rows = []
        removed_detail_rows = []

        for s_name, b_path in zip(samples_list, blast_files):
            short_s = re.sub(r"_junto$", "", s_name)
            sample_dir = os.path.dirname(os.path.dirname(b_path))
            cons_path = os.path.join(sample_dir, "Vsearch", f"{s_name}_consenso.fasta")

            cons_records = []
            if os.path.exists(cons_path):
                for rec in SeqIO.parse(cons_path, "fasta"):
                    qid = rec.id.split()[0]
                    sz_m = re.search(r"size=(\d+)", qid)
                    reads_v = int(sz_m.group(1)) if sz_m else 1
                    cons_records.append((qid, reads_v, len(rec.seq)))

            tot_cons_clusters = len(cons_records)
            tot_cons_reads = sum(r[1] for r in cons_records)
            no_blast_clusters, no_blast_reads = 0, 0
            non_target_clusters, non_target_reads = 0, 0
            trimal_drop_clusters, trimal_drop_reads = 0, 0
            retained_phylo_clusters, retained_phylo_reads = 0, 0

            no_blast_for_sample = []
            for qid, reads_v, raw_len in cons_records:
                uuid_m = re.search(r"centroid[_=]([0-9a-fA-F-]{8})", qid)
                short_u = uuid_m.group(1) if uuid_m else qid[:12]
                lca_rec = lca_records_by_key.get((s_name, qid))
                b_sp, b_acc, b_pident = blast_top_hit.get((s_name, qid), ("None", "N/A", 0.0))
                b_top_label = f"{b_sp} ({b_acc})" if b_acc and b_acc != "N/A" else b_sp

                if not lca_rec and (s_name, qid) not in blast_top_hit:
                    no_blast_clusters += 1
                    no_blast_reads += reads_v
                    seq_label = f"{short_s}_unmatched_{short_u}_size{reads_v}"
                    row_obj = {
                        "seq_id": seq_label,
                        "orig_qid": qid,
                        "sample": s_name,
                        "reads": reads_v,
                        "raw_len": raw_len,
                        "trimmed_12s_len": 0,
                        "stage": "1. Removed: No BLAST Hit (<pident% or Non-Mitochondrial)",
                        "blast_top_hit": "No Mitochondrial Hit",
                        "blast_pident": 0.0,
                        "lca_lineage": "Unclassified (No Hit in NCBI Mito RefSeq)",
                    }
                    all_audit_seq_rows.append(row_obj)
                    no_blast_for_sample.append(row_obj)
                elif (s_name, qid) not in qid_to_short_cluster:
                    non_target_clusters += 1
                    non_target_reads += reads_v
                    seq_label = f"{short_s}_nonTarget_{short_u}_size{reads_v}"
                    lin_str = lca_rec["lca_lineage"] if lca_rec else "Non-Target Lineage"
                    row_obj = {
                        "seq_id": seq_label,
                        "orig_qid": qid,
                        "sample": s_name,
                        "reads": reads_v,
                        "raw_len": raw_len,
                        "trimmed_12s_len": raw_len,
                        "stage": "2. Removed: Non-Target Clade (Outside Target TaxID in LCA)",
                        "blast_top_hit": b_top_label,
                        "blast_pident": b_pident,
                        "lca_lineage": lin_str,
                    }
                    all_audit_seq_rows.append(row_obj)
                    removed_detail_rows.append(row_obj)
                else:
                    fasta_id, c_short = qid_to_short_cluster[(s_name, qid)]
                    t_len = renamed_len_by_key.get((s_name, qid), raw_len)
                    lin_str = lca_rec["lca_lineage"] if lca_rec else "Target Clade"
                    phylo_rec = phylo_records_by_short.get((s_name, c_short))
                    if not phylo_rec:
                        trimal_drop_clusters += 1
                        trimal_drop_reads += reads_v
                        row_obj = {
                            "seq_id": fasta_id,
                            "orig_qid": qid,
                            "sample": s_name,
                            "reads": reads_v,
                            "raw_len": raw_len,
                            "trimmed_12s_len": t_len,
                            "stage": "3. Removed: trimAl / Alignment Coverage Filter (<60% non-gap coverage)",
                            "blast_top_hit": b_top_label,
                            "blast_pident": b_pident,
                            "lca_lineage": lin_str,
                        }
                        all_audit_seq_rows.append(row_obj)
                        removed_detail_rows.append(row_obj)
                    else:
                        retained_phylo_clusters += 1
                        retained_phylo_reads += reads_v
                        all_audit_seq_rows.append({
                            "seq_id": fasta_id,
                            "orig_qid": qid,
                            "sample": s_name,
                            "reads": reads_v,
                            "raw_len": raw_len,
                            "trimmed_12s_len": t_len,
                            "stage": f"4. Retained in Final Phylogeny (Closest: {phylo_rec['species']}, d={phylo_rec['p_dist']:.4f})",
                            "blast_top_hit": b_top_label,
                            "blast_pident": b_pident,
                            "lca_lineage": lin_str,
                        })

            # Include top 3 largest unmatched (No BLAST Hit) clusters per sample in the MultiQC detail table
            no_blast_for_sample.sort(key=lambda x: -x["reads"])
            for nb_row in no_blast_for_sample[:3]:
                removed_detail_rows.append(nb_row)

            blast_hit_clusters = tot_cons_clusters - no_blast_clusters
            blast_hit_reads = tot_cons_reads - no_blast_reads
            audit_summary_per_sample.append({
                "sample": s_name,
                "vsearch_clusters": tot_cons_clusters,
                "vsearch_reads": tot_cons_reads,
                "no_blast_clusters": no_blast_clusters,
                "no_blast_reads": no_blast_reads,
                "blast_hit_clusters": blast_hit_clusters,
                "blast_hit_reads": blast_hit_reads,
                "non_target_clusters": non_target_clusters,
                "non_target_reads": non_target_reads,
                "trimal_drop_clusters": trimal_drop_clusters,
                "trimal_drop_reads": trimal_drop_reads,
                "retained_clusters": retained_phylo_clusters,
                "retained_reads": retained_phylo_reads,
                "pct_phylo_reads": (retained_phylo_reads / tot_cons_reads * 100.0) if tot_cons_reads > 0 else 0.0,
            })

        if filtered_audit_tsv:
            os.makedirs(os.path.dirname(filtered_audit_tsv), exist_ok=True)
            with open(filtered_audit_tsv, "w") as fat:
                fat.write(
                    "Sample\tSequence_ID\tOriginal_VSEARCH_Header\tReads\tRaw_Consensus_Len_bp\tTrimmed_12S_Len_bp\t"
                    "Pipeline_Stage_Status\tBLAST_Top_Hit\tBLAST_Identity_Pct\tLCA_Lineage\n"
                )
                for r in all_audit_seq_rows:
                    fat.write(
                        f"{r['sample']}\t{r['seq_id']}\t{r['orig_qid']}\t{r['reads']}\t{r['raw_len']}\t{r['trimmed_12s_len']}\t"
                        f"{r['stage']}\t{r['blast_top_hit']}\t{r['blast_pident']:.2f}\t{r['lca_lineage']}\n"
                    )

        if filter_summary_mqc:
            os.makedirs(os.path.dirname(filter_summary_mqc), exist_ok=True)
            with open(filter_summary_mqc, "w") as fsm:
                fsm.write("# id: sequence_filtering_audit\n")
                fsm.write("# section_name: '1.2 Consensus Sequence Filtering & Retention Audit (VSEARCH -> BLAST -> Target Clade -> trimAl)'\n")
                fsm.write(
                    "# description: '<b>Methodology:</b> Tracks the fate of every non-chimeric VSEARCH consensus cluster across all downstream filtering checkpoints: "
                    "<b>(1) No BLAST Hit</b> (consensus sequences without any match in the NCBI Mitochondrial RefSeq database above the configured identity/coverage thresholds, representing off-target amplification or environmental background); "
                    "<b>(2) Removed: Non-Target Clade</b> (sequences that matched a mitochondrial reference in BLAST, such as plant mitochondria <i>Viridiplantae / Streptophyta</i>, but whose TaxonKit LCA lineage lies outside the configured phylogenetic target TaxID); "
                    "<b>(3) Removed: trimAl / Alignment Coverage</b> (target-clade fragments with &lt;60% non-gap coverage across the trimmed 12S alignment, excluded from tree reconstruction to preserve alignment quality and summarized in Section 1.3 and Section 4.1); and "
                    "<b>(4) Retained in Final Phylogeny</b> (high-confidence target-clade 12S consensus sequences with &ge;60% alignment coverage placed in the final FastTree ML tree).'\n"
                )
                fsm.write("# plot_type: 'table'\n")
                fsm.write("# pconfig:\n")
                fsm.write("#   id: 'sequence_filtering_audit_table'\n")
                fsm.write("#   title: 'Consensus Cluster & Read Retention Audit Across Pipeline Checkpoints'\n")
                fsm.write("#   col1_header: 'Sample'\n")
                fsm.write("# headers:\n")
                fsm.write("#   VSEARCH_Clusters:\n")
                fsm.write("#     title: 'VSEARCH Clusters'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   VSEARCH_Reads:\n")
                fsm.write("#     title: 'VSEARCH Reads'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Removed_No_BLAST_Clusters:\n")
                fsm.write("#     title: 'Removed: No BLAST Hit (Clusters)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Removed_No_BLAST_Reads:\n")
                fsm.write("#     title: 'Removed: No BLAST Hit (Reads)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   BLAST_Hit_Clusters:\n")
                fsm.write("#     title: 'Matched Mito DB (Clusters)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Removed_Non_Target_Clusters:\n")
                fsm.write("#     title: 'Removed: Non-Target Clade (Clusters)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Removed_Non_Target_Reads:\n")
                fsm.write("#     title: 'Removed: Non-Target Clade (Reads)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Removed_trimAl_Clusters:\n")
                fsm.write("#     title: 'Removed: trimAl (Clusters)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Removed_trimAl_Reads:\n")
                fsm.write("#     title: 'Removed: trimAl (Reads)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Retained_Phylo_Clusters:\n")
                fsm.write("#     title: 'Retained in Phylogeny (Clusters)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Retained_Phylo_Reads:\n")
                fsm.write("#     title: 'Retained in Phylogeny (Reads)'\n")
                fsm.write("#     format: '{:,.0f}'\n")
                fsm.write("#   Pct_Reads_in_Phylo:\n")
                fsm.write("#     title: '% Clustered Reads in Phylogeny'\n")
                fsm.write("#     suffix: '%'\n")
                fsm.write("#     format: '{:,.2f}'\n")
                fsm.write(
                    "Sample\tVSEARCH_Clusters\tVSEARCH_Reads\tRemoved_No_BLAST_Clusters\tRemoved_No_BLAST_Reads\t"
                    "BLAST_Hit_Clusters\tRemoved_Non_Target_Clusters\tRemoved_Non_Target_Reads\t"
                    "Removed_trimAl_Clusters\tRemoved_trimAl_Reads\tRetained_Phylo_Clusters\tRetained_Phylo_Reads\tPct_Reads_in_Phylo\n"
                )
                for a in audit_summary_per_sample:
                    fsm.write(
                        f"{a['sample']}\t{a['vsearch_clusters']}\t{a['vsearch_reads']}\t"
                        f"{a['no_blast_clusters']}\t{a['no_blast_reads']}\t{a['blast_hit_clusters']}\t"
                        f"{a['non_target_clusters']}\t{a['non_target_reads']}\t"
                        f"{a['trimal_drop_clusters']}\t{a['trimal_drop_reads']}\t"
                        f"{a['retained_clusters']}\t{a['retained_reads']}\t{a['pct_phylo_reads']:.2f}\n"
                    )

        if removed_seqs_mqc:
            os.makedirs(os.path.dirname(removed_seqs_mqc), exist_ok=True)
            removed_detail_rows.sort(key=lambda x: (x["stage"], x["sample"], -x["reads"]))
            with open(removed_seqs_mqc, "w") as rsm:
                rsm.write("# id: removed_sequences_detail\n")
                rsm.write("# section_name: '1.3 Detailed List of Removed / Filtered Sequences (Non-Target Clade, trimAl & Top Unmatched Clusters)'\n")
                rsm.write(
                    "# description: '<b>Methodology:</b> Lists individual consensus sequences removed during downstream filtering: "
                    "all sequences removed for belonging to a <b>Non-Target Clade</b> (e.g., plant mitochondrial hits outside <code>Vertebrata</code>), "
                    "all sequences removed by <b>trimAl / Alignment Coverage</b>, and the top 3 largest <b>No BLAST Hit</b> consensus clusters per sample. "
                    "The complete inventory of all 419 consensus sequences across all samples is saved in <code>all_samples/phylogeny/summaries/Filtered_and_Removed_Sequences_*.tsv</code>.'\n"
                )
                rsm.write("# plot_type: 'table'\n")
                rsm.write("# pconfig:\n")
                rsm.write("#   id: 'removed_sequences_detail_table'\n")
                rsm.write("#   title: 'Individual Sequences Removed by Non-Target Clade Filter, trimAl, or Lack of Mitochondrial BLAST Hit'\n")
                rsm.write("#   col1_header: 'Sequence_ID'\n")
                rsm.write("# headers:\n")
                rsm.write("#   Sample:\n")
                rsm.write("#     title: 'Sample'\n")
                rsm.write("#   Reads:\n")
                rsm.write("#     title: 'Reads'\n")
                rsm.write("#     format: '{:,.0f}'\n")
                rsm.write("#   Raw_Len_bp:\n")
                rsm.write("#     title: 'Consensus Len (bp)'\n")
                rsm.write("#     format: '{:,.0f}'\n")
                rsm.write("#   Trimmed_12S_Len_bp:\n")
                rsm.write("#     title: 'Aligned 12S Len (bp)'\n")
                rsm.write("#     format: '{:,.0f}'\n")
                rsm.write("#   Removal_Reason:\n")
                rsm.write("#     title: 'Removal Stage & Reason'\n")
                rsm.write("#   BLAST_Top_Hit:\n")
                rsm.write("#     title: 'BLAST Top Hit'\n")
                rsm.write("#   BLAST_Identity_Pct:\n")
                rsm.write("#     title: 'BLAST Identity (%)'\n")
                rsm.write("#     format: '{:,.2f}'\n")
                rsm.write("#   LCA_Lineage:\n")
                rsm.write("#     title: 'TaxonKit LCA Lineage'\n")
                rsm.write(
                    "Sequence_ID\tSample\tReads\tRaw_Len_bp\tTrimmed_12S_Len_bp\tRemoval_Reason\tBLAST_Top_Hit\tBLAST_Identity_Pct\tLCA_Lineage\n"
                )
                for r in removed_detail_rows:
                    rsm.write(
                        f"{r['seq_id']}\t{r['sample']}\t{r['reads']}\t{r['raw_len']}\t{r['trimmed_12s_len']}\t"
                        f"{r['stage']}\t{r['blast_top_hit']}\t{r['blast_pident']:.2f}\t{r['lca_lineage']}\n"
                    )

    def calc_shannon(counts_dict):
        tot = sum(counts_dict.values())
        if tot <= 0:
            return 0.0
        h = 0.0
        for c in counts_dict.values():
            if c > 0:
                p = c / float(tot)
                h -= p * math.log(p)
        return h

    samples_sorted = sorted(samples_list)
    all_lca_ranks = defaultdict(int)
    all_lca_genera = set()
    all_lca_species = set()
    all_phylo_ranks = defaultdict(int)
    all_phylo_species_counts = defaultdict(int)
    all_phylo_genera = set()
    all_phylo_species = set()

    for s in samples_sorted:
        for k, v in lca_sample_stats[s]["rank_counts"].items():
            all_lca_ranks[k] += v
        all_lca_genera |= lca_sample_stats[s]["genera"]
        all_lca_species |= lca_sample_stats[s]["species"]
        for k, v in phylo_sample_stats[s]["rank_counts"].items():
            all_phylo_ranks[k] += v
        for k, v in phylo_sample_stats[s]["species_counts"].items():
            all_phylo_species_counts[k] += v
        all_phylo_genera |= phylo_sample_stats[s]["genera"]
        all_phylo_species |= phylo_sample_stats[s]["species"]

    total_lca_reads = sum(lca_sample_stats[s]["reads"] for s in samples_sorted)
    total_lca_clusters = sum(lca_sample_stats[s]["clusters"] for s in samples_sorted)
    total_phylo_reads = sum(phylo_sample_stats[s]["reads"] for s in samples_sorted)
    total_phylo_clusters = sum(phylo_sample_stats[s]["clusters"] for s in samples_sorted)

    # 6. Write Section 4.2: Sample Taxonomic Richness & Shannon Diversity
    with open(richness_mqc, "w") as rf:
        rf.write("# id: sample_richness_comparison\n")
        rf.write("# section_name: '4.2 Taxonomic Richness & Diversity per Sample (BLAST+LCA vs. Phylogeny)'\n")
        rf.write(
            f"# description: '<b>Methodology:</b> Compares alpha-diversity and taxonomic richness metrics per sample/barcode (and pooled across the entire dataset in <code>ALL_SAMPLES_POOLED</code>) obtained via <b>BLAST + TaxonKit LCA</b> versus <b>FastTree Phylogenetic Placement</b>. "
            f"<b>Columns &amp; Interpretation:</b> <b>LCA_Reads / Phylo_Reads</b> and <b>LCA_Clusters / Phylo_Clusters</b> show total classified reads and consensus clusters; "
            f"<b>Richness_{rank_label}</b>, <b>Richness_Genus</b>, and <b>Richness_Species</b> report the exact number of unique taxa detected at the configured rank (<code>--lca-rank {defined_rank}</code>), Genus, and Species levels (note that LCA Species/Genus richness can be lower when equidistant BLAST hits collapse to higher taxonomic ranks, whereas Phylogeny always resolves to the closest reference tip); "
            f"<b>Shannon_H</b> is the Shannon-Wiener diversity index (<i>H&prime; = &minus;&sum; p<sub>i</sub> ln p<sub>i</sub></i>), which accounts for both taxon richness and read evenness.'\n"
        )
        rf.write("# plot_type: 'table'\n")
        rf.write("# pconfig:\n")
        rf.write("#   id: 'sample_richness_comparison_table'\n")
        rf.write(f"#   title: 'Sample Richness & Shannon Diversity: BLAST+LCA vs. Phylogeny (Defined Rank: {rank_label})'\n")
        rf.write("#   col1_header: 'Sample'\n")
        rf.write(
            f"Sample\tLCA_Reads\tLCA_Clusters\tLCA_Richness_{rank_label}\tLCA_Richness_Genus\tLCA_Richness_Species\tLCA_Shannon_H\t"
            f"Phylo_Reads\tPhylo_Clusters\tPhylo_Richness_{rank_label}\tPhylo_Richness_Genus\tPhylo_Richness_Species\tPhylo_Shannon_H\n"
        )
        for s in samples_sorted:
            lst = lca_sample_stats[s]
            pst = phylo_sample_stats[s]
            l_h = calc_shannon(lst["rank_counts"])
            p_h = calc_shannon(pst["species_counts"])
            rf.write(
                f"{s}\t{lst['reads']}\t{lst['clusters']}\t{len(lst['rank_counts'])}\t{len(lst['genera'])}\t{len(lst['species'])}\t{l_h:.3f}\t"
                f"{pst['reads']}\t{pst['clusters']}\t{len(pst['rank_counts'])}\t{len(pst['genera'])}\t{len(pst['species'])}\t{p_h:.3f}\n"
            )
        tot_l_h = calc_shannon(all_lca_ranks)
        tot_p_h = calc_shannon(all_phylo_species_counts)
        rf.write(
            f"ALL_SAMPLES_POOLED\t{total_lca_reads}\t{total_lca_clusters}\t{len(all_lca_ranks)}\t{len(all_lca_genera)}\t{len(all_lca_species)}\t{tot_l_h:.3f}\t"
            f"{total_phylo_reads}\t{total_phylo_clusters}\t{len(all_phylo_ranks)}\t{len(all_phylo_genera)}\t{len(all_phylo_species)}\t{tot_p_h:.3f}\n"
        )

    # 7. Write Section 4.3: Overall Taxonomic Abundance by Defined Rank
    all_rank_taxa = sorted(
        set(lca_rank_agg.keys()) | set(phylo_rank_agg.keys()),
        key=lambda t: -(lca_rank_agg[t]["reads"] + phylo_rank_agg[t]["reads"]),
    )

    with open(rank_abund_mqc, "w") as raf:
        raf.write("# id: overall_abundance_defined_rank\n")
        raf.write(f"# section_name: '4.3 Overall Taxonomic Abundance by {rank_label} (BLAST+LCA vs. Phylogeny)'\n")
        raf.write(
            f"# description: '<b>Methodology:</b> Dataset-wide synthesis comparing total read abundance, cluster counts, and sample occupancy (number of positive barcodes) for every <b>{rank_label}</b> (configured via <code>--lca-rank {defined_rank}</code>) identified by <b>BLAST+LCA</b> and by <b>Phylogenetic Placement</b> (mapped through the full NCBI lineage of each closest reference species). "
            f"<b>Columns &amp; Interpretation:</b> <b>Detected_By</b> indicates methodological concordance (<code>Both (LCA &amp; Phylogeny)</code> confirms high-confidence detection across both methods; <code>BLAST+LCA Only</code> occurs when a taxon falls outside the target phylogenetic clade or fails alignment coverage thresholds; <code>Phylogeny Only</code> occurs when LCA collapsed ambiguous BLAST hits to a higher rank above {rank_label}); "
            f"<b>Phylo_Species_Richness</b> reports the number of distinct closest reference species detected within that {rank_label} across the tree, and <b>Phylo_Mean_Distance</b> gives their read-weighted mean patristic genetic distance (substitutions/site).'\n"
        )
        raf.write("# plot_type: 'table'\n")
        raf.write("# pconfig:\n")
        raf.write("#   id: 'overall_abundance_defined_rank_table'\n")
        raf.write(f"#   title: 'Overall Abundance & Richness by {rank_label} (BLAST+LCA vs. Phylogeny)'\n")
        raf.write(f"#   col1_header: '{rank_label}'\n")
        raf.write(
            f"{rank_label}\tClass\tDetected_By\tLCA_Reads\tLCA_Relative_%\tLCA_Clusters\tLCA_Samples_Positive\t"
            f"Phylo_Reads\tPhylo_Relative_%\tPhylo_Clusters\tPhylo_Samples_Positive\tPhylo_Species_Richness\tPhylo_Mean_Distance\n"
        )
        for t_name in all_rank_taxa:
            la = lca_rank_agg.get(t_name, {"reads": 0, "clusters": 0, "samples": set(), "class": "Unclassified"})
            pa = phylo_rank_agg.get(t_name, {
                "reads": 0, "clusters": 0, "samples": set(), "species_set": set(),
                "dist_weighted_sum": 0.0, "class": "Unclassified"
            })
            l_reads = la["reads"]
            p_reads = pa["reads"]
            l_pct = (l_reads / float(total_lca_reads) * 100.0) if total_lca_reads > 0 else 0.0
            p_pct = (p_reads / float(total_phylo_reads) * 100.0) if total_phylo_reads > 0 else 0.0
            if l_reads > 0 and p_reads > 0:
                det_by = "Both (LCA & Phylogeny)"
            elif l_reads > 0:
                det_by = "BLAST+LCA Only"
            else:
                det_by = "Phylogeny Only"
            cls_name = pa["class"] if pa["class"] != "Unclassified" else la["class"]
            p_mean_d = (pa["dist_weighted_sum"] / float(p_reads)) if p_reads > 0 else 0.0
            clean_key = re.sub(r"[^A-Za-z0-9_.-]+", "_", t_name).strip("_") or "Unclassified"
            raf.write(
                f"{clean_key}\t{cls_name}\t{det_by}\t{l_reads}\t{l_pct:.2f}\t{la['clusters']}\t{len(la['samples'])}\t"
                f"{p_reads}\t{p_pct:.2f}\t{pa['clusters']}\t{len(pa['samples'])}\t{len(pa['species_set'])}\t{p_mean_d:.5f}\n"
            )

    # 8. Write Section 4.4: Overall Species-Level Abundance & Richness
    all_species_taxa = sorted(
        set(lca_species_agg.keys()) | set(phylo_species_agg.keys()),
        key=lambda t: -(lca_species_agg[t]["reads"] + phylo_species_agg[t]["reads"]),
    )

    with open(species_abund_mqc, "w") as saf, open(full_tsv, "w") as ftf:
        saf.write("# id: overall_species_abundance_comparison\n")
        saf.write("# section_name: '4.4 Overall Species-Level Abundance & Richness (BLAST+LCA vs. Phylogeny)'\n")
        saf.write(
            "# description: '<b>Methodology:</b> Complete dataset-wide species and terminal-taxon inventory across all samples, comparing direct species assignments from <b>BLAST+LCA</b> against closest evolutionary reference assignments from <b>FastTree Phylogeny</b>. "
            "<b>Columns &amp; Interpretation:</b> <b>Higher_Lineage</b> displays the taxonomic hierarchy (<code>Class &gt; Order &gt; Family</code>); "
            "<b>Detected_By</b> flags whether the species/taxon was identified by both approaches or uniquely by one; "
            "<b>LCA_Reads / Phylo_Reads</b>, <b>Relative_%</b>, <b>Clusters</b>, and <b>Samples_Positive</b> quantify abundance and prevalence across barcodes; "
            "<b>Phylo_Min_Distance</b> and <b>Phylo_Mean_Distance</b> (substitutions/site) evaluate species-level assignment confidence (values near <code>0.000&ndash;0.025</code> confirm exact or near-exact species matches, whereas higher values indicate an unrepresented species assigned to its closest congeneric/confamilial reference in NCBI or a co-amplified NUMT).'\n"
        )
        saf.write("# plot_type: 'table'\n")
        saf.write("# pconfig:\n")
        saf.write("#   id: 'overall_species_abundance_comparison_table'\n")
        saf.write("#   title: 'Overall Species Inventory & Abundance (BLAST+LCA vs. Phylogeny)'\n")
        saf.write("#   col1_header: 'Taxon_or_Species'\n")
        header_cols = (
            "Taxon_or_Species\tHigher_Lineage\tDetected_By\tLCA_Reads\tLCA_Relative_%\tLCA_Clusters\tLCA_Samples_Positive\t"
            "Phylo_Reads\tPhylo_Relative_%\tPhylo_Clusters\tPhylo_Samples_Positive\tPhylo_Min_Distance\tPhylo_Mean_Distance\n"
        )
        saf.write(header_cols)
        ftf.write(header_cols)

        for sp_name in all_species_taxa:
            lsa = lca_species_agg.get(sp_name, {"reads": 0, "clusters": 0, "samples": set(), "lineage": ""})
            psa = phylo_species_agg.get(sp_name, {
                "reads": 0, "clusters": 0, "samples": set(),
                "min_dist": 0.0, "dist_weighted_sum": 0.0, "lineage": ""
            })
            l_reads = lsa["reads"]
            p_reads = psa["reads"]
            l_pct = (l_reads / float(total_lca_reads) * 100.0) if total_lca_reads > 0 else 0.0
            p_pct = (p_reads / float(total_phylo_reads) * 100.0) if total_phylo_reads > 0 else 0.0
            if l_reads > 0 and p_reads > 0:
                det_by = "Both (LCA & Phylogeny)"
            elif l_reads > 0:
                det_by = "BLAST+LCA Only"
            else:
                det_by = "Phylogeny Only"
            lin_str = psa["lineage"] if psa["lineage"] else (lsa["lineage"] or "Unclassified")
            p_min_d = psa["min_dist"] if (p_reads > 0 and psa["min_dist"] != float("inf")) else 0.0
            p_mean_d = (psa["dist_weighted_sum"] / float(p_reads)) if p_reads > 0 else 0.0
            clean_sp = re.sub(r"[^A-Za-z0-9_.-]+", "_", sp_name).strip("_") or "Unclassified"
            row_line = (
                f"{clean_sp}\t{lin_str}\t{det_by}\t{l_reads}\t{l_pct:.2f}\t{lsa['clusters']}\t{len(lsa['samples'])}\t"
                f"{p_reads}\t{p_pct:.2f}\t{psa['clusters']}\t{len(psa['samples'])}\t{p_min_d:.5f}\t{p_mean_d:.5f}\n"
            )
            saf.write(row_line)
            ftf.write(row_line)
