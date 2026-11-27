############################################################################
#                    Phylogeny FASTA Preparation Helpers                   #
#                          MSc. Matheus Cosentino                          #
############################################################################

import os
import re
from Bio import SeqIO


def rename_sample_fasta_headers(lca_path, query_fasta, renamed_fasta, sample, phylo_target, config):
    """Filter sample consensus clusters matching target clade in LCA lineage and assign clean headers."""
    os.makedirs(os.path.dirname(renamed_fasta), exist_ok=True)

    target_id = int(phylo_target)
    phylo_names_map = {int(k): v for k, v in config.get("phylo_names", {}).items()}
    target_name = phylo_names_map.get(target_id, str(target_id))
    short_sample = re.sub(r"_junto$", "", sample)

    matched_info = {}
    if os.path.exists(lca_path):
        with open(lca_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 15:
                    continue
                query_id = parts[0]
                if query_id in matched_info:
                    continue
                full_lineage = parts[-2]
                reformatted_lineage = parts[-1]
                lineage_taxa = [t.strip() for t in (full_lineage + ";" + reformatted_lineage).split(";") if t.strip()]
                if target_name in lineage_taxa:
                    try:
                        qstart = int(parts[6])
                        qend = int(parts[7])
                        sstart = int(parts[8])
                        send = int(parts[9])
                    except ValueError:
                        qstart, qend, sstart, send = None, None, None, None
                    matched_info[query_id] = (qstart, qend, sstart, send)

    if not matched_info or not os.path.exists(query_fasta):
        with open(renamed_fasta, "w"):
            pass
        return

    records = []
    for record in SeqIO.parse(query_fasta, "fasta"):
        if record.id in matched_info:
            qstart, qend, sstart, send = matched_info[record.id]
            if qstart is not None and qend is not None and 1 <= qstart < qend <= len(record.seq):
                record.seq = record.seq[qstart - 1 : qend]
            if sstart is not None and send is not None and sstart > send:
                record.seq = record.seq.reverse_complement()
            records.append(record)

    def get_size(rec):
        m = re.search(r"size=(\d+)", rec.id + " " + rec.description)
        return int(m.group(1)) if m else 1

    records.sort(key=lambda r: (-get_size(r), r.id))

    with open(renamed_fasta, "w") as out_fh:
        for idx, record in enumerate(records, start=1):
            base_id = record.id.split()[0]
            size_val = get_size(record)
            original_desc = record.description if record.description else base_id
            record.id = f"{short_sample}_C{idx:02d}_size{size_val}"
            record.description = f"{sample}|C{idx:02d}|size={size_val}|{original_desc}"
            SeqIO.write(record, out_fh, "fasta")


def combine_renamed_fastas(input_fastas, combined_fasta):
    """Combine renamed sample FASTAs across all barcodes into a single non-redundant FASTA."""
    os.makedirs(os.path.dirname(combined_fasta), exist_ok=True)
    seen_ids = {}
    with open(combined_fasta, "w") as out_fh:
        for fasta_file in input_fastas:
            if os.path.exists(fasta_file) and os.path.getsize(fasta_file) > 0:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    base_id = record.id
                    if base_id in seen_ids:
                        seen_ids[base_id] += 1
                        record.id = f"{base_id}_{seen_ids[base_id]}"
                    else:
                        seen_ids[base_id] = 1
                    SeqIO.write(record, out_fh, "fasta")


def split_all_barcodes_fastas_by_group(combined_fasta, phylo_dist_files, lca_files, ref_tax_tsv, out_group_dir, rank):
    """
    Split combined consensus sequences across all samples into separate FASTA files
    by their identified taxonomic group at the configured taxonomic rank (e.g. family, order, class, genus, species),
    as well as subfolders by closest phylogenetic species and by BLAST+LCA rank.
    """
    import shutil

    if os.path.exists(out_group_dir):
        shutil.rmtree(out_group_dir)
    rank_key = str(rank).lower() if rank else "order"
    by_species_dir = os.path.join(out_group_dir, "by_species")
    by_lca_dir = os.path.join(out_group_dir, f"by_lca_{rank_key}")
    os.makedirs(by_species_dir, exist_ok=True)
    os.makedirs(by_lca_dir, exist_ok=True)

    species_to_ranks = {}
    genus_to_ranks = {}
    if ref_tax_tsv and os.path.exists(ref_tax_tsv):
        with open(ref_tax_tsv, "r") as tf:
            for line in tf:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 7:
                    _, p, c, o, f, g, s = parts[:7]
                    r_dict = {
                        "phylum": p if p and p != "Unclassified" else "Chordata",
                        "class": c if c and c != "Unclassified" else "",
                        "order": o if o and o != "Unclassified" else "",
                        "family": f if f and f != "Unclassified" else "",
                        "genus": g if g and g != "Unclassified" else "",
                        "species": s if s and s != "Unclassified" else "",
                    }
                    if s and s not in species_to_ranks:
                        species_to_ranks[s] = r_dict
                    if g and g not in genus_to_ranks:
                        genus_to_ranks[g] = r_dict

    # Parse BLAST+LCA lineage files to map (sample, original_qseqid) -> LCA rank dict
    lca_by_orig_qid = {}
    if lca_files:
        for lfile in lca_files:
            if not os.path.exists(lfile):
                continue
            with open(lfile, "r") as lf:
                for line in lf:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    qid = parts[0]
                    if qid in lca_by_orig_qid:
                        continue
                    full_lin = parts[-2] if len(parts) >= 2 else ""
                    ref_lin = parts[-1]
                    taxa = [t.strip() for t in ref_lin.split(";")]
                    while len(taxa) < 6:
                        taxa.append("")
                    if not taxa[1] and "Sauropsida" in full_lin:
                        taxa[1] = "Sauropsida"
                    lca_by_orig_qid[qid] = {
                        "phylum": taxa[0] or "Chordata",
                        "class": taxa[1] or "",
                        "order": taxa[2] or "",
                        "family": taxa[3] or "",
                        "genus": taxa[4] or "",
                        "species": taxa[5] or "",
                    }

    # Parse genetic_distance.tsv files:
    # Columns: Sample(0) | Cluster_ID(1) | Reads(2) | Closest_Species(3) | Closest_Accession(4) |
    #          Closest_Lineage(5) | Putative_Taxonomic_Level(6) | Patristic_Distance(7) | Terminal_Branch_Length(8)
    cluster_meta = {}
    for pfile in phylo_dist_files:
        if os.path.exists(pfile):
            with open(pfile, "r") as pf:
                header = pf.readline()
                for line in pf:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 8:
                        sample_name, cid, reads, c_sp, acc, c_lin, put_lvl, p_dist = parts[:8]
                    elif len(parts) >= 7:
                        sample_name, cid, reads, _, c_sp, acc, p_dist = parts[:7]
                        c_lin, put_lvl = "", ""
                    else:
                        continue
                    short_sample = re.sub(r"_junto$", "", sample_name)
                    meta_obj = {
                        "sample": sample_name,
                        "reads": reads,
                        "species": c_sp,
                        "lineage": c_lin,
                        "level": put_lvl,
                        "accession": acc,
                        "distance": p_dist,
                    }
                    # Index by full combined_fasta ID (e.g. barcode01_C01_size2518), sample_name_cid, and short cid
                    cluster_meta[f"{short_sample}_{cid}"] = meta_obj
                    cluster_meta[f"{sample_name}_{cid}"] = meta_obj
                    if cid.startswith(f"{short_sample}_") or cid.startswith(f"{sample_name}_"):
                        cluster_meta[cid] = meta_obj

    rank_handles = {}
    species_handles = {}
    lca_handles = {}

    def sanitize_filename(name):
        cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
        return cleaned.strip("_") or "Unclassified"

    def resolve_phylo_rank_taxon(c_sp, c_lin, lca_dict, r_key):
        if r_key == "species" and c_sp and c_sp != "Unclassified":
            return c_sp
        if r_key == "phylum":
            return "Chordata"
        gen_word = c_sp.split()[0] if c_sp and c_sp != "Unclassified" else ""
        sp_ranks = species_to_ranks.get(c_sp) or genus_to_ranks.get(gen_word) or {}
        val = sp_ranks.get(r_key, "")
        if val and val != "Unclassified":
            return val
        if c_lin and c_lin not in ("Unresolved", "Unclassified"):
            lin_parts = [x.strip() for x in c_lin.split(">")]
            # Closest_Lineage format: Class > Order > Family > Genus
            idx_map = {"class": 0, "order": 1, "family": 2, "genus": 3}
            if r_key in idx_map and len(lin_parts) > idx_map[r_key] and lin_parts[idx_map[r_key]]:
                return lin_parts[idx_map[r_key]]
        # Fallback to BLAST+LCA rank if sequence was not placed in tree
        if lca_dict:
            l_val = lca_dict.get(r_key, "")
            if l_val:
                return l_val
            for higher_r in ("family", "order", "class", "phylum"):
                if lca_dict.get(higher_r):
                    return f"Unclassified_{lca_dict[higher_r]}"
        return "Unclassified"

    def resolve_lca_rank_taxon(lca_dict, r_key):
        if not lca_dict:
            return "Unclassified"
        val = lca_dict.get(r_key, "")
        if val:
            return val
        for higher_r in ("genus", "family", "order", "class", "phylum"):
            if lca_dict.get(higher_r):
                return f"Unclassified_{lca_dict[higher_r]}"
        return "Unclassified"

    try:
        if os.path.exists(combined_fasta) and os.path.getsize(combined_fasta) > 0:
            for record in SeqIO.parse(combined_fasta, "fasta"):
                cid = record.id.split()[0]
                desc_parts = record.description.split("|", 3)
                orig_qid = desc_parts[3].split()[0] if len(desc_parts) >= 4 else ""
                lca_dict = lca_by_orig_qid.get(orig_qid, {})

                meta = cluster_meta.get(cid, {})
                c_sp = meta.get("species", "") or lca_dict.get("species", "") or "Unclassified"
                c_lin = meta.get("lineage", "")
                put_lvl = meta.get("level", "BLAST+LCA")
                acc = meta.get("accession", "N/A")
                p_dist = meta.get("distance", "N/A")

                rank_taxon = resolve_phylo_rank_taxon(c_sp, c_lin, lca_dict, rank_key)
                lca_taxon = resolve_lca_rank_taxon(lca_dict, rank_key)

                record.description = (
                    f"{cid} | phylo_{rank_key}={rank_taxon} | lca_{rank_key}={lca_taxon} "
                    f"| closest_species={c_sp} ({acc}) | distance={p_dist} | level={put_lvl}"
                )

                r_fname = sanitize_filename(rank_taxon)
                s_fname = sanitize_filename(c_sp)
                l_fname = sanitize_filename(lca_taxon)

                if r_fname not in rank_handles:
                    rank_handles[r_fname] = open(os.path.join(out_group_dir, f"{r_fname}.fasta"), "w")
                SeqIO.write(record, rank_handles[r_fname], "fasta")

                if s_fname not in species_handles:
                    species_handles[s_fname] = open(os.path.join(by_species_dir, f"{s_fname}.fasta"), "w")
                SeqIO.write(record, species_handles[s_fname], "fasta")

                if l_fname not in lca_handles:
                    lca_handles[l_fname] = open(os.path.join(by_lca_dir, f"{l_fname}.fasta"), "w")
                SeqIO.write(record, lca_handles[l_fname], "fasta")
    finally:
        for h in rank_handles.values():
            h.close()
        for h in species_handles.values():
            h.close()
        for h in lca_handles.values():
            h.close()


