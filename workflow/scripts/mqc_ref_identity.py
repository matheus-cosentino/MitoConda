############################################################################
#       Reference DB Pairwise Identity & Genetic Distance Histogram        #
#                          MSc. Matheus Cosentino                          #
############################################################################

import os
import re
import json
import numpy as np
from collections import defaultdict
from Bio import SeqIO


def run_ref_identity_histogram_mqc(
    ref_aln, tax_tsv, dist_tsv, summary_mqc, hist_svg, mqc_html,
    phylo_target, gene, config, rank=None
):
    """Compute vectorized pairwise identity & genetic distance on the aligned reference DB and build MultiQC visualizations."""
    os.makedirs(os.path.dirname(dist_tsv), exist_ok=True)

    if rank:
        defined_rank = str(rank).lower()
    else:
        cfg_ranks = config.get("taxonkit", {}).get("lca_ranks") or config.get("taxonkit", {}).get("lca_rank", ["order"])
        if isinstance(cfg_ranks, str):
            cfg_ranks = [cfg_ranks]
        defined_rank = str(cfg_ranks[0]).lower() if cfg_ranks else "order"
    if defined_rank not in ("species", "genus", "family", "order", "class", "phylum"):
        defined_rank = "order"

    cfg_pidents = config.get("blast", {}).get("perc_identity", [90])
    if not isinstance(cfg_pidents, list):
        cfg_pidents = [cfg_pidents]

    tax_map = {}
    if os.path.exists(tax_tsv):
        with open(tax_tsv, "r") as tf:
            for line in tf:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 7:
                    sid, p, c, o, f, g, s = parts[:7]
                    tax_map[sid] = {
                        "phylum": p if p and p != "Unclassified" else "Chordata",
                        "class": c if c and c != "Unclassified" else "",
                        "order": o if o and o != "Unclassified" else "",
                        "family": f if f and f != "Unclassified" else "",
                        "genus": g if g and g != "Unclassified" else "",
                        "species": s if s and s != "Unclassified" else "",
                    }

    records = list(SeqIO.parse(ref_aln, "fasta")) if os.path.exists(ref_aln) else []
    n_seqs = len(records)

    if n_seqs < 2:
        for p in (dist_tsv, summary_mqc, hist_svg, mqc_html):
            with open(p, "w") as fh:
                fh.write("")
        return

    aln_len = max(len(r.seq) for r in records)
    seq_ids = [r.id.split()[0] for r in records]

    arr = np.full((n_seqs, aln_len), ord("-"), dtype=np.uint8)
    for idx, r in enumerate(records):
        s_bytes = np.frombuffer(str(r.seq).upper().encode("ascii", errors="replace"), dtype=np.uint8)
        arr[idx, : len(s_bytes)] = s_bytes[:aln_len]

    valid_col_mask = (
        (arr == ord("A")) | (arr == ord("C")) | (arr == ord("G")) | (arr == ord("T"))
    ).sum(axis=0) >= 2
    arr = arr[:, valid_col_mask]

    is_A = (arr == ord("A")).astype(np.float32)
    is_C = (arr == ord("C")).astype(np.float32)
    is_G = (arr == ord("G")).astype(np.float32)
    is_T = (arr == ord("T")).astype(np.float32)
    valid_bases = is_A + is_C + is_G + is_T

    shared_valid = valid_bases @ valid_bases.T
    exact_matches = (is_A @ is_A.T) + (is_C @ is_C.T) + (is_G @ is_G.T) + (is_T @ is_T.T)
    del is_A, is_C, is_G, is_T, valid_bases, arr

    min_overlap = 200.0
    with np.errstate(divide="ignore", invalid="ignore"):
        pident_mat = np.where(
            shared_valid >= min_overlap,
            (exact_matches / np.maximum(shared_valid, 1.0)) * 100.0,
            np.nan,
        )
    np.fill_diagonal(pident_mat, np.nan)
    del shared_valid, exact_matches

    sp_arr = np.array([tax_map.get(sid, {}).get("species", "") for sid in seq_ids], dtype=object)
    gen_arr = np.array([tax_map.get(sid, {}).get("genus", "") for sid in seq_ids], dtype=object)
    fam_arr = np.array([tax_map.get(sid, {}).get("family", "") for sid in seq_ids], dtype=object)
    ord_arr = np.array([tax_map.get(sid, {}).get("order", "") for sid in seq_ids], dtype=object)
    cls_arr = np.array([tax_map.get(sid, {}).get("class", "") for sid in seq_ids], dtype=object)
    phy_arr = np.array([tax_map.get(sid, {}).get("phylum", "Chordata") for sid in seq_ids], dtype=object)

    iu, ju = np.triu_indices(n_seqs, k=1)
    pair_vals = pident_mat[iu, ju]
    valid_pair_mask = ~np.isnan(pair_vals)
    iu = iu[valid_pair_mask]
    ju = ju[valid_pair_mask]
    pair_vals = pair_vals[valid_pair_mask]

    same_sp = (sp_arr[iu] != "") & (sp_arr[iu] == sp_arr[ju])
    same_gen = (gen_arr[iu] != "") & (gen_arr[iu] == gen_arr[ju])
    same_fam = (fam_arr[iu] != "") & (fam_arr[iu] == fam_arr[ju])
    same_ord = (ord_arr[iu] != "") & (ord_arr[iu] == ord_arr[ju])
    same_cls = (cls_arr[iu] != "") & (cls_arr[iu] == cls_arr[ju])
    same_phy = (phy_arr[iu] != "") & (phy_arr[iu] == phy_arr[ju])

    hier_masks = {
        "Species (Intra-Species)": same_sp,
        "Genus (Same Genus, Diff. Species)": same_gen & (~same_sp),
        "Family (Same Family, Diff. Genus)": same_fam & (~same_gen),
        "Order (Same Order, Diff. Family)": same_ord & (~same_fam),
        "Class (Same Class, Diff. Order)": same_cls & (~same_ord),
        "Phylum (Same Phylum, Diff. Class)": ~same_cls,
    }
    hier_data = {k: pair_vals[m] for k, m in hier_masks.items() if np.any(m)}

    nn_lists = defaultdict(list)
    for i in range(n_seqs):
        row = pident_mat[i]
        v_mask = ~np.isnan(row)
        if not np.any(v_mask):
            continue
        s_i, g_i, f_i, o_i, c_i = sp_arr[i], gen_arr[i], fam_arr[i], ord_arr[i], cls_arr[i]

        m_sp = v_mask & (sp_arr == s_i) if s_i else np.zeros(n_seqs, dtype=bool)
        m_gen = v_mask & (gen_arr == g_i) & (sp_arr != s_i) if g_i else np.zeros(n_seqs, dtype=bool)
        m_fam = v_mask & (fam_arr == f_i) & (gen_arr != g_i) if f_i else np.zeros(n_seqs, dtype=bool)
        m_ord = v_mask & (ord_arr == o_i) & (fam_arr != f_i) if o_i else np.zeros(n_seqs, dtype=bool)
        m_cls = v_mask & (cls_arr == c_i) & (ord_arr != o_i) if c_i else np.zeros(n_seqs, dtype=bool)
        m_phy = v_mask & (cls_arr != c_i) if c_i else np.zeros(n_seqs, dtype=bool)

        if np.any(m_sp):
            nn_lists["Species (Best Hit Same Species)"].append(float(np.nanmax(row[m_sp])))
        if np.any(m_gen):
            nn_lists["Genus (Best Hit Same Genus)"].append(float(np.nanmax(row[m_gen])))
        if np.any(m_fam):
            nn_lists["Family (Best Hit Same Family)"].append(float(np.nanmax(row[m_fam])))
        if np.any(m_ord):
            nn_lists["Order (Best Hit Same Order)"].append(float(np.nanmax(row[m_ord])))
        if np.any(m_cls):
            nn_lists["Class (Best Hit Same Class)"].append(float(np.nanmax(row[m_cls])))
        if np.any(m_phy):
            nn_lists["Phylum (Best Hit Diff. Class)"].append(float(np.nanmax(row[m_phy])))

    nn_data = {k: np.array(v, dtype=np.float32) for k, v in nn_lists.items() if len(v) > 0}

    if "Species (Intra-Species)" not in hier_data and "Genus (Best Hit Same Genus)" in nn_data:
        ordered_hier = {"Species (Nearest Sister Species in Genus)": nn_data["Genus (Best Hit Same Genus)"]}
        ordered_hier.update(hier_data)
        hier_data = ordered_hier

    # Build whole-database distributions for each taxonomic rank (Phylum, Class, Order, Family, Genus, Species)
    # so 100% of the reference database is represented (Global All-Pairs, Global Exclusive Rank,
    # Global Best-Hit Nearest Neighbor, Inter-Rank Barcoding Gap, and Partitioned by Vertebrate Class),
    # plus full per-taxon statistics across all taxa in the database (without arbitrary top-12 truncation).
    rank_meta = {
        "phylum": {
            "arr": phy_arr,
            "same_mask": same_phy,
            "excl_mask": ~same_cls,
            "excl_label": "Entire DB: Same Phylum, Diff. Class",
            "nn_key": "Phylum (Best Hit Diff. Class)",
            "gap_mask": same_cls,
            "gap_label": "Entire DB: Intra-Class Pairs (Within-Class Baseline)",
        },
        "class": {
            "arr": cls_arr,
            "same_mask": same_cls,
            "excl_mask": same_cls & (~same_ord),
            "excl_label": "Entire DB: Same Class, Diff. Order",
            "nn_key": "Class (Best Hit Same Class)",
            "gap_mask": ~same_cls,
            "gap_label": "Entire DB: Inter-Class Pairs (Diff. Class)",
        },
        "order": {
            "arr": ord_arr,
            "same_mask": same_ord,
            "excl_mask": same_ord & (~same_fam),
            "excl_label": "Entire DB: Same Order, Diff. Family",
            "nn_key": "Order (Best Hit Same Order)",
            "gap_mask": same_cls & (~same_ord),
            "gap_label": "Entire DB: Inter-Order Pairs (Same Class, Diff. Order)",
        },
        "family": {
            "arr": fam_arr,
            "same_mask": same_fam,
            "excl_mask": same_fam & (~same_gen),
            "excl_label": "Entire DB: Same Family, Diff. Genus",
            "nn_key": "Family (Best Hit Same Family)",
            "gap_mask": same_ord & (~same_fam),
            "gap_label": "Entire DB: Inter-Family Pairs (Same Order, Diff. Family)",
        },
        "genus": {
            "arr": gen_arr,
            "same_mask": same_gen,
            "excl_mask": same_gen & (~same_sp),
            "excl_label": "Entire DB: Same Genus, Diff. Species",
            "nn_key": "Genus (Best Hit Same Genus)",
            "gap_mask": same_fam & (~same_gen),
            "gap_label": "Entire DB: Inter-Genus Pairs (Same Family, Diff. Genus)",
        },
    }

    rank_group_data = {}
    all_taxa_per_rank = {}

    for r_name, r_info in rank_meta.items():
        r_arr = r_info["arr"]
        r_same_mask = r_info["same_mask"]
        r_title = r_name.title()
        r_dict = {}

        # 1. Whole-database aggregate series for this rank
        if np.any(r_same_mask):
            r_dict[f"Entire DB: All Intra-{r_title} Pairs"] = pair_vals[r_same_mask]
        if np.any(r_info["excl_mask"]):
            r_dict[r_info["excl_label"]] = pair_vals[r_info["excl_mask"]]
        if r_info["nn_key"] in nn_data:
            r_dict[f"Entire DB: Best Hit / Nearest Neighbor in {r_title}"] = nn_data[r_info["nn_key"]]
        if np.any(r_info["gap_mask"]):
            r_dict[r_info["gap_label"]] = pair_vals[r_info["gap_mask"]]

        # 2. Partition all Intra-Rank pairs across every Vertebrate Class in the database (covers 100% of DB)
        sub_iu = iu[r_same_mask]
        sub_vals = pair_vals[r_same_mask]
        if r_name != "phylum" and len(sub_vals) > 0:
            sub_classes = cls_arr[sub_iu]
            cls_buckets = defaultdict(list)
            for c_val, p_val in zip(sub_classes, sub_vals):
                if c_val:
                    cls_buckets[c_val].append(float(p_val))
            for c_name, c_vals in sorted(cls_buckets.items(), key=lambda kv: (-len(kv[1]), kv[0])):
                if len(c_vals) >= 1:
                    r_dict[f"Class {c_name}: All Intra-{r_title} Pairs"] = np.array(c_vals, dtype=np.float32)

        # 3. Full per-taxon statistics across ALL taxa of this rank in the entire database (no top-12 cutoff)
        sub_taxa = r_arr[sub_iu]
        group_buckets = defaultdict(list)
        for t_val, p_val in zip(sub_taxa, sub_vals):
            if t_val:
                group_buckets[t_val].append(float(p_val))
        sorted_groups = sorted(group_buckets.items(), key=lambda kv: (-len(kv[1]), kv[0]))
        all_taxa_per_rank[r_name] = [
            (g_name, np.array(g_vals, dtype=np.float32))
            for g_name, g_vals in sorted_groups
            if len(g_vals) >= 1
        ]

        rank_group_data[r_name] = r_dict

    sp_tab_dict = {}
    if np.any(same_sp):
        sp_tab_dict["Entire DB: Intra-Species Pairs"] = pair_vals[same_sp]
    if "Genus (Best Hit Same Genus)" in nn_data:
        sp_tab_dict["Entire DB: Nearest Sister Species in Genus (Best Hit)"] = nn_data["Genus (Best Hit Same Genus)"]
    if np.any(same_gen & (~same_sp)):
        sp_tab_dict["Entire DB: All Congeneric Sister-Species Pairs"] = pair_vals[same_gen & (~same_sp)]
    for i in range(n_seqs):
        row = pident_mat[i]
        v_mask = ~np.isnan(row)
        s_i, g_i, c_i = sp_arr[i], gen_arr[i], cls_arr[i]
        if g_i and c_i:
            m_gen = v_mask & (gen_arr == g_i) & (sp_arr != s_i)
            if np.any(m_gen):
                sp_tab_dict.setdefault(f"Class {c_i}: Nearest Sister Species", []).append(float(np.nanmax(row[m_gen])))
    for k in list(sp_tab_dict.keys()):
        if isinstance(sp_tab_dict[k], list):
            if len(sp_tab_dict[k]) >= 1:
                sp_tab_dict[k] = np.array(sp_tab_dict[k], dtype=np.float32)
            else:
                del sp_tab_dict[k]
    rank_group_data["species"] = sp_tab_dict

    id_bin_edges = np.linspace(45.0, 100.0, 111)
    id_bin_centers = ((id_bin_edges[:-1] + id_bin_edges[1:]) / 2.0).round(3).tolist()

    dist_bin_edges = np.linspace(0.0, 0.55, 111)
    dist_bin_centers = ((dist_bin_edges[:-1] + dist_bin_edges[1:]) / 2.0).round(4).tolist()

    def summarize_array(vals):
        n_v = int(len(vals))
        if n_v == 0:
            return None
        d_vals = np.clip(1.0 - (vals / 100.0), 0.0, 1.0)
        id_counts, _ = np.histogram(vals, bins=id_bin_edges)
        id_pcts = (id_counts / float(n_v) * 100.0).round(3).tolist()
        dist_counts, _ = np.histogram(d_vals, bins=dist_bin_edges)
        dist_pcts = (dist_counts / float(n_v) * 100.0).round(3).tolist()
        return {
            "n": n_v,
            "mean": float(np.mean(vals)),
            "sd": float(np.std(vals)),
            "median": float(np.median(vals)),
            "p05": float(np.percentile(vals, 5)),
            "p25": float(np.percentile(vals, 25)),
            "p75": float(np.percentile(vals, 75)),
            "p95": float(np.percentile(vals, 95)),
            "min": float(np.min(vals)),
            "max": float(np.max(vals)),
            "d_mean": float(np.mean(d_vals)),
            "d_sd": float(np.std(d_vals)),
            "d_median": float(np.median(d_vals)),
            "d_p05": float(np.percentile(d_vals, 5)),
            "d_p25": float(np.percentile(d_vals, 25)),
            "d_p75": float(np.percentile(d_vals, 75)),
            "d_p95": float(np.percentile(d_vals, 95)),
            "d_min": float(np.min(d_vals)),
            "d_max": float(np.max(d_vals)),
            "counts": id_counts.astype(int).tolist(),
            "pcts": id_pcts,
            "d_counts": dist_counts.astype(int).tolist(),
            "d_pcts": dist_pcts,
        }

    summary_rows = []
    with open(dist_tsv, "w") as dfh:
        dfh.write(
            "View_Category\tTaxon_or_Rank\tN_Comparisons\tMean_Identity\tSD_Identity\t"
            "Median_Identity\tP05_Identity\tP25_Identity\tP75_Identity\tP95_Identity\tMin_Identity\tMax_Identity\t"
            "Mean_Genetic_Distance\tMedian_Genetic_Distance\tP25_Genetic_Distance\tP75_Genetic_Distance\n"
        )
        for label, arr_v in hier_data.items():
            st = summarize_array(arr_v)
            if st:
                dfh.write(
                    f"Hierarchy_All_Pairs\t{label}\t{st['n']}\t{st['mean']:.2f}\t{st['sd']:.2f}\t"
                    f"{st['median']:.2f}\t{st['p05']:.2f}\t{st['p25']:.2f}\t{st['p75']:.2f}\t"
                    f"{st['p95']:.2f}\t{st['min']:.2f}\t{st['max']:.2f}\t"
                    f"{st['d_mean']:.4f}\t{st['d_median']:.4f}\t{st['d_p25']:.4f}\t{st['d_p75']:.4f}\n"
                )
                summary_rows.append(("All_Pairs_" + label, "Hierarchy (All Pairs)", label, st))

        for label, arr_v in nn_data.items():
            st = summarize_array(arr_v)
            if st:
                dfh.write(
                    f"Hierarchy_Best_Hit\t{label}\t{st['n']}\t{st['mean']:.2f}\t{st['sd']:.2f}\t"
                    f"{st['median']:.2f}\t{st['p05']:.2f}\t{st['p25']:.2f}\t{st['p75']:.2f}\t"
                    f"{st['p95']:.2f}\t{st['min']:.2f}\t{st['max']:.2f}\t"
                    f"{st['d_mean']:.4f}\t{st['d_median']:.4f}\t{st['d_p25']:.4f}\t{st['d_p75']:.4f}\n"
                )
                summary_rows.append(("Best_Hit_" + label, "Best Hit (Nearest Neighbor)", label, st))

        for r_key in ("phylum", "class", "order", "family", "genus", "species"):
            for g_label, arr_v in rank_group_data.get(r_key, {}).items():
                st = summarize_array(arr_v)
                if st:
                    dfh.write(
                        f"Whole_DB_{r_key.title()}\t{g_label}\t{st['n']}\t{st['mean']:.2f}\t{st['sd']:.2f}\t"
                        f"{st['median']:.2f}\t{st['p05']:.2f}\t{st['p25']:.2f}\t{st['p75']:.2f}\t"
                        f"{st['p95']:.2f}\t{st['min']:.2f}\t{st['max']:.2f}\t"
                        f"{st['d_mean']:.4f}\t{st['d_median']:.4f}\t{st['d_p25']:.4f}\t{st['d_p75']:.4f}\n"
                    )
                    summary_rows.append((f"WholeDB_{r_key}_{g_label}", f"Whole DB — {r_key.title()}", g_label, st))

            # Also write every individual taxon in this rank across the entire database to dist_tsv (and to 5.2 for Class/Order/Family)
            for t_name, arr_v in all_taxa_per_rank.get(r_key, []):
                st = summarize_array(arr_v)
                if st:
                    dfh.write(
                        f"Taxon_{r_key.title()}\t{t_name}\t{st['n']}\t{st['mean']:.2f}\t{st['sd']:.2f}\t"
                        f"{st['median']:.2f}\t{st['p05']:.2f}\t{st['p25']:.2f}\t{st['p75']:.2f}\t"
                        f"{st['p95']:.2f}\t{st['min']:.2f}\t{st['max']:.2f}\t"
                        f"{st['d_mean']:.4f}\t{st['d_median']:.4f}\t{st['d_p25']:.4f}\t{st['d_p75']:.4f}\n"
                    )
                    if r_key in ("class", "order", "family") or r_key == defined_rank:
                        summary_rows.append((f"Taxon_{r_key}_{t_name}", f"Intra-{r_key.title()} (All DB Taxa)", t_name, st))

    with open(summary_mqc, "w") as sfh:
        sfh.write("# id: ref_db_identity_summary\n")
        sfh.write(f"# section_name: '5.2 Reference DB Genetic Distance & Identity Summary (Phylum, Class, Order, Family, Genus, Species)'\n")
        sfh.write(
            "# description: '<b>Methodology:</b> Summarizes the pairwise genetic distance (<i>p</i>-distance = <code>1 &minus; Identity/100</code>, substitutions/site) and sequence identity (%) computed across all aligned reference pairs (&ge;200 bp non-gap overlap) in <code>"
            f"{phylo_target}_{gene}_ref_aln.fasta</code> (N={n_seqs} reference sequences), stratified across <b>Phylum</b>, <b>Class</b>, <b>Order</b>, <b>Family</b>, <b>Genus</b>, and <b>Species</b>. "
            "<b>Categories &amp; Interpretation:</b> <b>Hierarchy (All Pairs)</b> reports statistics across all pairwise combinations sharing that lowest taxonomic rank; "
            "<b>Best Hit (Nearest Neighbor)</b> reports the distribution of the single closest relative for each reference sequence at that rank (modeling expected BLAST Top Hit distances when a species is absent from NCBI); "
            "<b>Intra-Phylum / Intra-Class / Intra-Order</b> break down genetic distances within major vertebrate classes and orders.'\n"
        )
        sfh.write("# plot_type: 'table'\n")
        sfh.write("# pconfig:\n")
        sfh.write("#    id: 'ref_db_identity_summary_table'\n")
        sfh.write("#    title: 'Reference Database Genetic Distance & Sequence Identity (%) by Taxonomic Level'\n")
        sfh.write("#    col1_header: 'Comparison_ID'\n")
        sfh.write("# headers:\n")
        sfh.write("#    Category:\n")
        sfh.write("#      title: 'Category'\n")
        sfh.write("#    Taxonomic_Level_or_Group:\n")
        sfh.write("#      title: 'Taxonomic Level / Lineage Group'\n")
        sfh.write("#    N_Comparisons:\n")
        sfh.write("#      title: 'N Pairs'\n")
        sfh.write("#      format: '{:,.0f}'\n")
        sfh.write("#    Mean_Genetic_Dist:\n")
        sfh.write("#      title: 'Mean Genetic Distance'\n")
        sfh.write("#      format: '{:,.4f}'\n")
        sfh.write("#    Median_Genetic_Dist:\n")
        sfh.write("#      title: 'Median Genetic Distance'\n")
        sfh.write("#      format: '{:,.4f}'\n")
        sfh.write("#    IQR_Genetic_Dist:\n")
        sfh.write("#      title: 'IQR Genetic Distance (P25 - P75)'\n")
        sfh.write("#    Mean_Identity_Pct:\n")
        sfh.write("#      title: 'Mean Identity (%)'\n")
        sfh.write("#      format: '{:,.2f}'\n")
        sfh.write("#    Median_Identity_Pct:\n")
        sfh.write("#      title: 'Median Identity (%)'\n")
        sfh.write("#      format: '{:,.2f}'\n")
        sfh.write("#    IQR_Identity_Pct:\n")
        sfh.write("#      title: 'IQR Identity % (P25 - P75)'\n")
        sfh.write(
            "Comparison_ID\tCategory\tTaxonomic_Level_or_Group\tN_Comparisons\tMean_Genetic_Dist\tMedian_Genetic_Dist\tIQR_Genetic_Dist\t"
            "Mean_Identity_Pct\tMedian_Identity_Pct\tIQR_Identity_Pct\n"
        )
        for row_id, cat_name, grp_name, st in summary_rows:
            clean_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", row_id).strip("_")
            sfh.write(
                f"{clean_id}\t{cat_name}\t{grp_name}\t{st['n']}\t"
                f"{st['d_mean']:.4f}\t{st['d_median']:.4f}\t{st['d_p25']:.4f} - {st['d_p75']:.4f}\t"
                f"{st['mean']:.2f}\t{st['median']:.2f}\t{st['p25']:.2f} - {st['p75']:.2f}\n"
            )

    palette = [
        "#e11d48", "#2563eb", "#16a34a", "#d97706", "#7c3aed",
        "#0891b2", "#db2777", "#4f46e5", "#059669", "#ea580c", "#64748b",
    ]
    svg_w, svg_h = 1100, 760
    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_w} {svg_h}" width="100%" height="100%" style="background:#ffffff; font-family:Arial,sans-serif;">',
        f'<text x="25" y="32" font-size="18" font-weight="bold" fill="#0f172a">Reference Database 12S Pairwise Identity Distribution ({n_seqs} Reference Sequences)</text>',
        f'<text x="25" y="52" font-size="12" fill="#475569">Top: Comparison across Taxonomic Levels | Bottom: Distribution within Defined Taxonomic Rank ({defined_rank.title()})</text>',
    ]

    def render_svg_panel(data_dict, top_y, panel_h, title_str):
        left_x, right_x = 75, svg_w - 280
        pw = right_x - left_x
        ph = panel_h - 65
        py = top_y + 30
        svg_lines.append(f'<text x="{left_x}" y="{top_y + 18}" font-size="14" font-weight="bold" fill="#1e293b">{title_str}</text>')
        svg_lines.append(f'<rect x="{left_x}" y="{py}" width="{pw}" height="{ph}" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1"/>')

        all_stats = [(k, summarize_array(v)) for k, v in data_dict.items()]
        all_stats = [(k, s) for k, s in all_stats if s is not None]
        max_pct = max((max(s["pcts"]) for _, s in all_stats), default=10.0)
        max_pct = max(max_pct * 1.12, 5.0)

        for x_val in range(45, 101, 5):
            px = left_x + ((x_val - 45.0) / 55.0) * pw
            svg_lines.append(f'<line x1="{px:.1f}" y1="{py}" x2="{px:.1f}" y2="{py + ph}" stroke="#e2e8f0" stroke-width="1"/>')
            svg_lines.append(f'<text x="{px:.1f}" y="{py + ph + 16}" font-size="10" fill="#475569" text-anchor="middle">{x_val}%</text>')

        for step_i in range(5):
            y_val = (max_pct / 4.0) * step_i
            py_tick = py + ph - (y_val / max_pct) * ph
            svg_lines.append(f'<line x1="{left_x}" y1="{py_tick:.1f}" x2="{right_x}" y2="{py_tick:.1f}" stroke="#e2e8f0" stroke-width="1"/>')
            svg_lines.append(f'<text x="{left_x - 8}" y="{py_tick + 3:.1f}" font-size="10" fill="#475569" text-anchor="end">{y_val:.1f}%</text>')

        for idx, (label, st) in enumerate(all_stats):
            color = palette[idx % len(palette)]
            pts = [f"{left_x:.1f},{py + ph:.1f}"]
            for b_center, pct_v in zip(id_bin_centers, st["pcts"]):
                bx = left_x + ((b_center - 45.0) / 55.0) * pw
                by = py + ph - (pct_v / max_pct) * ph
                pts.append(f"{bx:.1f},{by:.1f}")
            pts.append(f"{right_x:.1f},{py + ph:.1f}")
            pts_str = " ".join(pts)
            svg_lines.append(f'<polygon points="{pts_str}" fill="{color}" fill-opacity="0.24" stroke="{color}" stroke-width="1.8"/>')

            ly = py + 16 + idx * 20
            if ly < py + ph + 15:
                svg_lines.append(f'<rect x="{right_x + 14}" y="{ly - 9}" width="12" height="12" rx="2" fill="{color}" fill-opacity="0.7" stroke="{color}"/>')
                short_lbl = label if len(label) <= 30 else label[:28] + ".."
                svg_lines.append(
                    f'<text x="{right_x + 32}" y="{ly + 1}" font-size="10.5" fill="#1e293b">{short_lbl} (med {st["median"]:.1f}%)</text>'
                )

        for pid_th in cfg_pidents:
            try:
                pth = float(pid_th)
                if 45.0 <= pth <= 100.0:
                    tx = left_x + ((pth - 45.0) / 55.0) * pw
                    svg_lines.append(f'<line x1="{tx:.1f}" y1="{py}" x2="{tx:.1f}" y2="{py + ph}" stroke="#dc2626" stroke-width="1.6" stroke-dasharray="5,4"/>')
                    svg_lines.append(f'<text x="{tx - 4:.1f}" y="{py + 14}" font-size="10" font-weight="bold" fill="#dc2626" text-anchor="end">BLAST {pth:.0f}%</text>')
            except ValueError:
                pass

    render_svg_panel(hier_data, 68, 330, "1. Pairwise Identity (%) Across Taxonomic Hierarchy Levels")
    render_svg_panel(
        rank_group_data.get(defined_rank, {}),
        415,
        330,
        f"2. Pairwise Identity (%) Within Defined Taxonomic Rank: {defined_rank.title()}",
    )
    svg_lines.append("</svg>")
    with open(hist_svg, "w") as svf:
        svf.write("\n".join(svg_lines))

    views_payload = {}
    view_definitions = [
        ("hier_all", "Overall Hierarchy (All Pairs)", hier_data),
        ("hier_nn", "Overall Hierarchy (Best Hit / Top Hit)", nn_data),
        ("phylum", "By Phylum (Chordata)", rank_group_data.get("phylum", {})),
        ("class", "By Class", rank_group_data.get("class", {})),
        ("order", "By Order", rank_group_data.get("order", {})),
        ("family", "By Family", rank_group_data.get("family", {})),
        ("genus", "By Genus", rank_group_data.get("genus", {})),
        ("species", "By Species", rank_group_data.get("species", {})),
    ]

    for v_key, v_title, v_dict in view_definitions:
        series_list = []
        for idx, (s_label, s_arr) in enumerate(v_dict.items()):
            st = summarize_array(s_arr)
            if not st:
                continue
            color = palette[idx % len(palette)]
            series_list.append({
                "name": f"{s_label} (med: {st['median']:.1f}%, d={st['d_median']:.3f}, n={st['n']:,})",
                "short_name": s_label,
                "color": color,
                "pcts": st["pcts"],
                "counts": st["counts"],
                "d_pcts": st["d_pcts"],
                "d_counts": st["d_counts"],
                "median": round(st["median"], 2),
                "mean": round(st["mean"], 2),
                "p25": round(st["p25"], 2),
                "p75": round(st["p75"], 2),
                "d_median": round(st["d_median"], 4),
                "d_p25": round(st["d_p25"], 4),
                "d_p75": round(st["d_p75"], 4),
                "n": st["n"],
            })
        views_payload[v_key] = {
            "title": v_title,
            "series": series_list,
        }

    div_id = f"mqc_ref_identity_hist_{phylo_target}_{gene}"
    payload_json = json.dumps({
        "id_bins": id_bin_centers,
        "dist_bins": dist_bin_centers,
        "defined_rank": defined_rank,
        "blast_pidents": [float(p) for p in cfg_pidents],
        "n_seqs": n_seqs,
        "views": views_payload,
    })

    html_content = f"""<!--
id: 'ref_db_identity_histogram'
section_name: '5.1 Reference DB Genetic Distance & Identity Distribution (Rank: {defined_rank.title()})'
description: '<b>Methodology:</b> Evaluates the taxonomic resolution and barcoding gap of the target marker across the aligned reference database (<code>{phylo_target}_{gene}_ref_aln.fasta</code>, N={n_seqs} reference sequences). Pairwise sequence identity (%) and genetic distance (<i>p</i>-distance = <code>1 &minus; Identity/100</code>, substitutions/site) are computed across all sequence pairs sharing &ge;200 bp of non-gap alignment and stratified by their NCBI taxonomic lineage. <b>Top Panel (Taxonomic Hierarchy Comparison):</b> Compares distributions across mutually exclusive shared ranks (<i>Species</i>, <i>Genus</i>, <i>Family</i>, <i>Order</i>, <i>Class</i>, and <i>Phylum</i>); toggle between <b>Hierarchy (All Pairs)</b> and <b>Hierarchy (Best Hit / BLAST Top Hit)</b> (which records each sequence&rsquo;s single closest relative at each rank, modeling expected BLAST Top Hit identities when a sample species is absent from NCBI). <b>Bottom Panel (Intra-Group Distribution):</b> Opens on the configured rank (<code>--lca-rank {defined_rank}</code>) and lets you switch interactively across <b>Phylum</b>, <b>Class</b>, <b>Order</b>, <b>Family</b>, <b>Genus</b>, and <b>Species</b> to inspect within-clade sequence divergence. Use the <b>Metric</b> buttons to switch the X-axis between <i>Identity (%)</i> and <i>Genetic Distance</i>, and the <b>Y-Axis</b> buttons to switch between normalized <i>Relative % per Level</i> and raw <i>Pair Count</i>. The red dashed line marks the active BLAST identity cutoff.'
plot_type: 'html'
-->
<style>
.mqc-custom-content-plot {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }}
.mqc-custom-content-plot .btn {{ cursor: pointer; border-radius: 4px; padding: 3px 8px; font-size: 12px; line-height: 1.4; border: 1px solid #cbd5e1; background: #f8fafc; color: #334155; transition: all 0.15s; }}
.mqc-custom-content-plot .btn:hover {{ background: #e2e8f0; }}
.mqc-custom-content-plot .btn-primary {{ background: #2563eb !important; color: #ffffff !important; border-color: #1d4ed8 !important; }}
</style>
<div class="mqc-custom-content-plot" style="border:1px solid #e2e8f0; border-radius:8px; padding:14px; background:#ffffff;">
    <div style="display:flex; flex-wrap:wrap; gap:8px; align-items:center; justify-content:space-between; margin-bottom:12px; padding-bottom:10px; border-bottom:1px solid #f1f5f9;">
        <div style="display:flex; flex-wrap:wrap; gap:6px; align-items:center;">
            <span style="font-size:12px; font-weight:700; color:#334155; margin-right:4px;">Taxonomic Level:</span>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="hier_all" style="font-weight:600;">Hierarchy (All Pairs)</button>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="hier_nn" style="font-weight:600;">Hierarchy (Best Hit / BLAST Top Hit)</button>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="phylum">Phylum</button>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="class">Class</button>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="order">Order</button>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="family">Family</button>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="genus">Genus</button>
            <button type="button" class="btn btn-xs btn-default ref-id-view-btn_{div_id}" data-view="species">Species</button>
        </div>
        <div style="display:flex; flex-wrap:wrap; gap:10px; align-items:center;">
            <div style="display:flex; gap:6px; align-items:center;">
                <span style="font-size:12px; font-weight:700; color:#334155;">Metric:</span>
                <button type="button" class="btn btn-xs btn-primary ref-id-metric-btn_{div_id}" data-metric="identity">Identity (%)</button>
                <button type="button" class="btn btn-xs btn-default ref-id-metric-btn_{div_id}" data-metric="distance">Genetic Distance</button>
            </div>
            <div style="display:flex; gap:6px; align-items:center;">
                <span style="font-size:12px; font-weight:700; color:#334155;">Y-Axis:</span>
                <button type="button" class="btn btn-xs btn-primary ref-id-norm-btn_{div_id}" data-mode="pct">Relative % per Level</button>
                <button type="button" class="btn btn-xs btn-default ref-id-norm-btn_{div_id}" data-mode="count">Pair Count</button>
            </div>
        </div>
    </div>
    <div id="{div_id}_top" style="width:100%; height:440px;"></div>
    <div id="{div_id}_bottom" style="width:100%; height:440px; margin-top:14px; border-top:1px dashed #e2e8f0; padding-top:12px;"></div>
</div>
<script type="text/javascript">
(function initMitoCondaRefIdentityHist() {{
    if (typeof Plotly === "undefined") {{
        if (!window._mitoconda_loading_plotly) {{
            window._mitoconda_loading_plotly = true;
            var s = document.createElement("script");
            s.src = "https://cdn.plot.ly/plotly-2.35.2.min.js";
            document.head.appendChild(s);
        }}
        setTimeout(initMitoCondaRefIdentityHist, 100);
        return;
    }}
    var payload = {payload_json};
    var topDiv = document.getElementById("{div_id}_top");
    var bottomDiv = document.getElementById("{div_id}_bottom");
    if (!topDiv || !bottomDiv) return;

    var currentTopView = "hier_all";
    var currentBottomView = payload.defined_rank || "order";
    var currentNormMode = "pct";
    var currentMetric = "identity";

    function buildPlotForView(targetEl, viewKey, panelPrefix) {{
        var viewObj = payload.views[viewKey] || payload.views["hier_all"];
        var isId = (currentMetric === "identity");
        var bins = isId ? payload.id_bins : payload.dist_bins;
        var isPct = (currentNormMode === "pct");
        var traces = [];

        for (var i = 0; i < viewObj.series.length; i++) {{
            var s = viewObj.series[i];
            var pctArr = isId ? s.pcts : s.d_pcts;
            var cntArr = isId ? s.counts : s.d_counts;
            var yVals = isPct ? pctArr : cntArr;
            var custom = [];
            for (var j = 0; j < bins.length; j++) {{
                custom.push([
                    pctArr[j],
                    cntArr[j],
                    s.median,
                    s.p25,
                    s.p75,
                    s.d_median,
                    s.d_p25,
                    s.d_p75,
                    s.n
                ]);
            }}
            traces.push({{
                type: "bar",
                name: s.name,
                x: bins,
                y: yVals,
                width: isId ? 0.48 : 0.0048,
                marker: {{
                    color: s.color,
                    opacity: 0.55,
                    line: {{ color: s.color, width: 1.1 }}
                }},
                customdata: custom,
                hovertemplate:
                    "<b>" + s.short_name + "</b><br>" +
                    (isId ? "Sequence Identity Bin: <b>%{{x:.2f}}%</b><br>" : "Genetic Distance Bin: <b>%{{x:.4f}} subs/site</b><br>") +
                    "Proportion in Level: <b>%{{customdata[0]:.2f}}%</b><br>" +
                    "Comparisons in Bin: <b>%{{customdata[1]:,}}</b><br>" +
                    "Group Median Identity: <b>%{{customdata[2]:.2f}}%</b> (IQR: %{{customdata[3]:.1f}}%–%{{customdata[4]:.1f}}%)<br>" +
                    "Group Median Distance: <b>%{{customdata[5]:.4f}}</b> (IQR: %{{customdata[6]:.4f}}–%{{customdata[7]:.4f}}, N=%{{customdata[8]:,}})" +
                    "<extra></extra>"
            }});
        }}

        var shapes = [];
        var annotations = [];
        for (var k = 0; k < payload.blast_pidents.length; k++) {{
            var thPid = payload.blast_pidents[k];
            var thVal = isId ? thPid : Math.max(0.0, 1.0 - (thPid / 100.0));
            shapes.push({{
                type: "line",
                x0: thVal,
                x1: thVal,
                y0: 0,
                y1: 1,
                yref: "paper",
                line: {{ color: "#dc2626", width: 2, dash: "dash" }}
            }});
            annotations.push({{
                x: thVal,
                y: 0.98,
                xref: "x",
                yref: "paper",
                text: isId ? ("BLAST Cutoff (" + thPid + "%)") : ("BLAST Cutoff (d=" + thVal.toFixed(2) + ")"),
                showarrow: false,
                xanchor: isId ? "right" : "left",
                font: {{ size: 11, color: "#dc2626" }},
                bgcolor: "rgba(255,255,255,0.85)"
            }});
        }}

        var metricLabel = isId ? "Sequence Identity (%)" : "Genetic Distance (substitutions/site)";
        var layout = {{
            title: {{
                text: panelPrefix + viewObj.title + " — " + metricLabel + " (" + (isPct ? "Relative Distribution %" : "Pair Counts") + ")",
                font: {{ size: 14, color: "#0f172a" }},
                x: 0.01
            }},
            barmode: "overlay",
            bargap: 0.02,
            dragmode: "zoom",
            hovermode: "closest",
            margin: {{ l: 65, r: 20, t: 45, b: 50 }},
            plot_bgcolor: "#ffffff",
            paper_bgcolor: "#ffffff",
            xaxis: {{
                title: metricLabel,
                range: isId ? [50, 100.5] : [-0.005, 0.50],
                dtick: isId ? 5 : 0.05,
                gridcolor: "#f1f5f9",
                zeroline: false
            }},
            yaxis: {{
                title: isPct ? "% of Comparisons in Taxonomic Level" : "Number of Comparisons (Pairs)",
                gridcolor: "#f1f5f9",
                zeroline: true
            }},
            shapes: shapes,
            annotations: annotations,
            legend: {{
                orientation: "v",
                x: 1.01,
                y: 1,
                bgcolor: "rgba(255,255,255,0.92)",
                bordercolor: "#e2e8f0",
                borderwidth: 1,
                font: {{ size: 11 }}
            }}
        }};

        var config = {{
            responsive: true,
            displaylogo: false,
            scrollZoom: true
        }};

        Plotly.react(targetEl, traces, layout, config);
    }}

    function updateButtonStyles() {{
        var viewBtns = document.querySelectorAll(".ref-id-view-btn_{div_id}");
        for (var i = 0; i < viewBtns.length; i++) {{
            var v = viewBtns[i].getAttribute("data-view");
            var isActive = (v === currentTopView || v === currentBottomView);
            viewBtns[i].className = "btn btn-xs " + (isActive ? "btn-primary" : "btn-default") + " ref-id-view-btn_{div_id}";
        }}
        var normBtns = document.querySelectorAll(".ref-id-norm-btn_{div_id}");
        for (var j = 0; j < normBtns.length; j++) {{
            var m = normBtns[j].getAttribute("data-mode");
            normBtns[j].className = "btn btn-xs " + (m === currentNormMode ? "btn-primary" : "btn-default") + " ref-id-norm-btn_{div_id}";
        }}
        var metricBtns = document.querySelectorAll(".ref-id-metric-btn_{div_id}");
        for (var k = 0; k < metricBtns.length; k++) {{
            var mt = metricBtns[k].getAttribute("data-metric");
            metricBtns[k].className = "btn btn-xs " + (mt === currentMetric ? "btn-primary" : "btn-default") + " ref-id-metric-btn_{div_id}";
        }}
    }}

    function renderAll() {{
        buildPlotForView(topDiv, currentTopView, "1. Taxonomic Level Comparison: ");
        buildPlotForView(bottomDiv, currentBottomView, "2. Intra-Group Distribution — ");
        updateButtonStyles();
    }}

    var viewBtns = document.querySelectorAll(".ref-id-view-btn_{div_id}");
    for (var i = 0; i < viewBtns.length; i++) {{
        viewBtns[i].addEventListener("click", function() {{
            var clickedView = this.getAttribute("data-view");
            if (clickedView === "hier_all" || clickedView === "hier_nn") {{
                currentTopView = clickedView;
            }} else {{
                currentBottomView = clickedView;
            }}
            renderAll();
        }});
    }}

    var normBtns = document.querySelectorAll(".ref-id-norm-btn_{div_id}");
    for (var j = 0; j < normBtns.length; j++) {{
        normBtns[j].addEventListener("click", function() {{
            currentNormMode = this.getAttribute("data-mode");
            renderAll();
        }});
    }}

    var metricBtns = document.querySelectorAll(".ref-id-metric-btn_{div_id}");
    for (var k = 0; k < metricBtns.length; k++) {{
        metricBtns[k].addEventListener("click", function() {{
            currentMetric = this.getAttribute("data-metric");
            renderAll();
        }});
    }}

    renderAll();
}})();
</script>
"""
    with open(mqc_html, "w") as mfh:
        mfh.write(html_content)
