############################################################################
#                       LCA Summary & FASTA Splitting                      #
#                          MSc. Matheus Cosentino                          #
############################################################################

import os
import re
import glob
from collections import defaultdict
from Bio import SeqIO


def summarize_lca_abundance(blast_out, abundance_table, rank):
    """Summarize read abundance per taxon at the requested taxonomic rank."""
    os.makedirs(os.path.dirname(abundance_table), exist_ok=True)

    rank_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
    rank_index = rank_map.get(rank, 4)

    taxon_abundance = defaultdict(int)
    processed_centroids = set()

    if os.path.exists(blast_out):
        with open(blast_out, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue

                centroid_id = parts[0]
                if centroid_id in processed_centroids:
                    continue
                processed_centroids.add(centroid_id)

                size_match = re.search(r"size=(\d+)", centroid_id)
                count = int(size_match.group(1)) if size_match else 1

                full_lineage = parts[-2] if len(parts) >= 2 else ""
                lineage = parts[-1]
                taxa = [t.strip() for t in lineage.split(";")]
                if len(taxa) > 1 and not taxa[1] and "Sauropsida" in full_lineage:
                    taxa[1] = "Sauropsida"
                non_empty_taxa = [t for t in taxa if t]

                if len(taxa) > rank_index and taxa[rank_index]:
                    taxon = taxa[rank_index]
                elif non_empty_taxa:
                    taxon = non_empty_taxa[-1]
                else:
                    taxon = "Unclassified"

                taxon_abundance[taxon] += count

    with open(abundance_table, "w") as out:
        out.write("Taxon\tAbundance\n")
        for taxon, total in sorted(taxon_abundance.items(), key=lambda x: x[1], reverse=True):
            out.write(f"{taxon}\t{total}\n")


def split_fastas_by_lca(lca_path, query_fasta, fasta_dir, rank):
    """Trim consensus sequences to aligned 12S region, orient strand, and split by LCA taxon."""
    rank_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
    rank_index = rank_map.get(rank, 4)

    seq_info = {}
    if os.path.exists(lca_path):
        with open(lca_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue

                query_id = parts[0]
                if query_id in seq_info:
                    continue

                try:
                    qstart = int(parts[6])
                    qend = int(parts[7])
                    sstart = int(parts[8])
                    send = int(parts[9])
                except ValueError:
                    qstart, qend, sstart, send = None, None, None, None

                full_lineage = parts[-2] if len(parts) >= 2 else ""
                lineage = parts[-1]
                taxa = [t.strip() for t in lineage.split(";")]
                if len(taxa) > 1 and not taxa[1] and "Sauropsida" in full_lineage:
                    taxa[1] = "Sauropsida"
                non_empty_taxa = [t for t in taxa if t]
                if len(taxa) > rank_index and taxa[rank_index]:
                    taxon = taxa[rank_index]
                elif non_empty_taxa:
                    taxon = non_empty_taxa[-1]
                else:
                    taxon = "Unclassified"

                taxon_clean = taxon.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
                seq_info[query_id] = (taxon_clean, qstart, qend, sstart, send)

    os.makedirs(fasta_dir, exist_ok=True)
    for old_fasta in glob.glob(os.path.join(fasta_dir, "*.fasta")):
        os.remove(old_fasta)

    handles = {}
    try:
        if os.path.exists(query_fasta):
            for record in SeqIO.parse(query_fasta, "fasta"):
                if record.id in seq_info:
                    t_name, qstart, qend, sstart, send = seq_info[record.id]
                    if qstart is not None and qend is not None and 1 <= qstart < qend <= len(record.seq):
                        record.seq = record.seq[qstart - 1 : qend]
                    if sstart is not None and send is not None and sstart > send:
                        record.seq = record.seq.reverse_complement()

                    file_path = os.path.join(fasta_dir, f"{t_name}.fasta")
                    if t_name not in handles:
                        handles[t_name] = open(file_path, "w")
                    SeqIO.write(record, handles[t_name], "fasta")
    finally:
        for h in handles.values():
            h.close()

