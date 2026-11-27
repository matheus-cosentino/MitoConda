############################################################################
#          Interactive Midpoint-Rooted Phylogenetic Tree (MultiQC)         #
#                          MSc. Matheus Cosentino                          #
############################################################################

import os
import re
import math
import json
import heapq
import html
from collections import defaultdict, deque
from Bio import Phylo


def run_phylo_tree_plot_mqc(tree_path, svg_file, mqc_html, samples, phylo_target, pident, gene):
    """Midpoint-root full integrated tree, prune to <=5 nodes from sample tips, and render SVG + Plotly HTML."""
    os.makedirs(os.path.dirname(svg_file), exist_ok=True)
    tree = Phylo.read(tree_path, "newick")

    sample_prefixes = tuple(
        p for s in samples for p in (f"{s}_", f"{re.sub(r'_junto$', '', s)}_")
    )

    cluster_tips = []
    ref_tips = []
    for tip in tree.get_terminals():
        tname = tip.name or ""
        if sample_prefixes and tname.startswith(sample_prefixes):
            cluster_tips.append(tip)
        else:
            ref_tips.append(tip)

    def midpoint_root_fast(phylo_tree):
        terminals = phylo_tree.get_terminals()
        if len(terminals) < 2:
            return phylo_tree

        parent_map = {}
        t_adj = defaultdict(list)
        c_by_id = {}
        rogue_clade_ids = set()
        for p in phylo_tree.find_clades(order="level"):
            c_by_id[id(p)] = p
            for ch in p.clades:
                c_by_id[id(ch)] = ch
                parent_map[id(ch)] = id(p)
                w = float(ch.branch_length) if ch.branch_length is not None else 0.0
                if w < 0.0:
                    w = 0.0
                if w > 0.50 or id(p) in rogue_clade_ids:
                    rogue_clade_ids.add(id(ch))
                t_adj[id(p)].append((id(ch), w))
                t_adj[id(ch)].append((id(p), w))

        clean_refs = [t for t in ref_tips if id(t) not in rogue_clade_ids]
        if len(clean_refs) >= 2:
            valid_term_ids = set(id(t) for t in clean_refs)
            start_leaf_id = id(clean_refs[0])
        elif len(ref_tips) >= 2:
            valid_term_ids = set(id(t) for t in ref_tips)
            start_leaf_id = id(ref_tips[0])
        else:
            valid_term_ids = set(id(t) for t in terminals)
            start_leaf_id = id(terminals[0])

        def bfs_farthest(start_node_id):
            dist = {start_node_id: 0.0}
            prev = {start_node_id: None}
            stack = [start_node_id]
            farthest_t = start_node_id
            max_d = -1.0
            while stack:
                u = stack.pop()
                du = dist[u]
                if u in valid_term_ids and du > max_d:
                    max_d = du
                    farthest_t = u
                for v, w in t_adj[u]:
                    if v not in dist:
                        dist[v] = du + w
                        prev[v] = u
                        stack.append(v)
            return farthest_t, max_d, dist, prev

        u_id, _, _, _ = bfs_farthest(start_leaf_id)
        v_id, max_diam, dist_from_u, prev_from_u = bfs_farthest(u_id)

        if max_diam <= 0.0:
            return phylo_tree

        half_d = max_diam / 2.0
        path = []
        curr = v_id
        while curr is not None:
            path.append(curr)
            curr = prev_from_u.get(curr)
        path.reverse()

        for i in range(len(path) - 1):
            a_id, b_id = path[i], path[i + 1]
            da, db = dist_from_u[a_id], dist_from_u[b_id]
            if da <= half_d <= db:
                edge_len = db - da
                if parent_map.get(a_id) == b_id:
                    child_clade = c_by_id[a_id]
                    outgroup_bl = max(0.0, min(edge_len, half_d - da))
                elif parent_map.get(b_id) == a_id:
                    child_clade = c_by_id[b_id]
                    dist_to_v_from_b = max_diam - db
                    outgroup_bl = max(0.0, min(edge_len, half_d - dist_to_v_from_b))
                else:
                    break
                phylo_tree.root_with_outgroup(child_clade, outgroup_branch_length=outgroup_bl)
                break
        return phylo_tree

    tree = midpoint_root_fast(tree)
    tree.root.branch_length = 0.0

    adj = defaultdict(list)
    u_adj = defaultdict(list)
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            w = float(child.branch_length) if child.branch_length is not None else 0.0
            if w < 0.0:
                w = 0.0
            adj[id(parent)].append((id(child), w))
            adj[id(child)].append((id(parent), w))
            u_adj[id(parent)].append(id(child))
            u_adj[id(child)].append(id(parent))

    nearest_ref_for_cluster = {}
    nearest_ref_id_for_cluster = {}
    if cluster_tips and ref_tips:
        dist_to_ref = {}
        nearest_ref_name = {}
        nearest_ref_id = {}
        pq = []
        for rtip in ref_tips:
            rid = id(rtip)
            rname = rtip.name or "Reference"
            dist_to_ref[rid] = 0.0
            nearest_ref_name[rid] = rname
            nearest_ref_id[rid] = rid
            heapq.heappush(pq, (0.0, rname, rid, rid))
        while pq:
            d, rname, r_origin_id, u = heapq.heappop(pq)
            if d > dist_to_ref.get(u, float("inf")):
                continue
            for v, w in adj[u]:
                nd = d + w
                if nd < dist_to_ref.get(v, float("inf")) or (
                    nd == dist_to_ref.get(v, float("inf")) and rname < nearest_ref_name.get(v, "")
                ):
                    dist_to_ref[v] = nd
                    nearest_ref_name[v] = rname
                    nearest_ref_id[v] = r_origin_id
                    heapq.heappush(pq, (nd, rname, r_origin_id, v))
        for ctip in cluster_tips:
            cid = id(ctip)
            rname = nearest_ref_name.get(cid, "Reference")
            m_r = re.match(r"^(.+)_([A-Za-z]{2}_?\d+\.\d+)(?:_\d+)?$", rname)
            clean_r = f"{m_r.group(1).replace('_', ' ')} ({m_r.group(2)})" if m_r else rname.replace("_", " ")
            nearest_ref_for_cluster[cid] = (clean_r, dist_to_ref.get(cid, 0.0))
            if cid in nearest_ref_id:
                nearest_ref_id_for_cluster[cid] = nearest_ref_id[cid]

    max_ref_node_dist = 5
    cluster_tip_ids = set(id(t) for t in cluster_tips)
    ref_tip_ids = set(id(t) for t in ref_tips)
    all_tip_ids = cluster_tip_ids | ref_tip_ids

    if cluster_tips:
        sample_only_clades = set()
        def mark_sample_only(clade):
            if not clade.clades:
                is_s = id(clade) in cluster_tip_ids
                if is_s:
                    sample_only_clades.add(id(clade))
                return is_s
            child_flags = [mark_sample_only(c) for c in clade.clades]
            if all(child_flags):
                sample_only_clades.add(id(clade))
                return True
            return False

        mark_sample_only(tree.root)

        node_dist = {}
        bfs_q = deque()
        for sid in sample_only_clades:
            node_dist[sid] = 0
            bfs_q.append(sid)

        keep_leaf_ids = set(cluster_tip_ids)
        while bfs_q:
            u = bfs_q.popleft()
            du = node_dist[u]
            for v in u_adj[u]:
                if v in ref_tip_ids:
                    if du <= max_ref_node_dist:
                        keep_leaf_ids.add(v)
                elif v not in all_tip_ids:
                    nd = 0 if v in sample_only_clades else (du + 1)
                    if nd <= max_ref_node_dist and nd < node_dist.get(v, 999999):
                        node_dist[v] = nd
                        bfs_q.append(v)

        for ctip in cluster_tips:
            n_rid = nearest_ref_id_for_cluster.get(id(ctip))
            if n_rid is not None:
                keep_leaf_ids.add(n_rid)

        def prune_clade_fast(clade, is_root=False):
            if not clade.clades:
                return clade if id(clade) in keep_leaf_ids else None
            kept_children = []
            for ch in clade.clades:
                res = prune_clade_fast(ch, is_root=False)
                if res is not None:
                    kept_children.append(res)
            if not kept_children:
                return None
            if len(kept_children) == 1 and not is_root:
                only_child = kept_children[0]
                p_bl = float(clade.branch_length) if clade.branch_length is not None else 0.0
                c_bl = float(only_child.branch_length) if only_child.branch_length is not None else 0.0
                only_child.branch_length = max(0.0, p_bl) + max(0.0, c_bl)
                return only_child
            clade.clades = kept_children
            return clade

        prune_clade_fast(tree.root, is_root=True)

    pruned_root = tree.root
    pruned_root.branch_length = 0.0
    pruned_tree = tree
    pruned_tree.ladderize()

    x_coords = {}
    y_coords = {}
    leaves = pruned_tree.get_terminals()
    n_leaves = max(1, len(leaves))

    for idx, leaf in enumerate(leaves):
        y_coords[id(leaf)] = idx

    def assign_x(clade, curr_x):
        bl = float(clade.branch_length) if clade.branch_length is not None else 0.0
        if bl < 0.0:
            bl = 0.0
        nx = curr_x + bl
        x_coords[id(clade)] = nx
        for child in clade.clades:
            assign_x(child, nx)

    def assign_y(clade):
        if not clade.clades:
            return y_coords[id(clade)]
        child_ys = [assign_y(c) for c in clade.clades]
        y_coords[id(clade)] = sum(child_ys) / len(child_ys)
        return y_coords[id(clade)]

    assign_x(pruned_root, 0.0)
    assign_y(pruned_root)

    max_x = max(x_coords.values()) if x_coords else 1.0
    if max_x <= 0.0:
        max_x = 1.0

    palette = [
        "#e11d48", "#2563eb", "#16a34a", "#d97706", "#9333ea",
        "#0891b2", "#ea580c", "#4f46e5", "#0d9488", "#c026d3",
    ]
    sample_colors = {}
    for idx, s in enumerate(sorted(samples)):
        short_s = re.sub(r"_junto$", "", s)
        color = palette[idx % len(palette)]
        sample_colors[s] = color
        sample_colors[short_s] = color

    def get_tip_style(tip_name):
        for s_key, col in sample_colors.items():
            if tip_name.startswith(f"{s_key}_"):
                m_sz = re.search(r"size[_=]?(\d+)", tip_name)
                reads = int(m_sz.group(1)) if m_sz else 1
                radius = min(7.5, max(3.5, 2.2 + math.log10(max(1, reads)) * 1.6))
                return True, col, radius, tip_name
        m_new = re.match(r"^(.+)_([A-Za-z]{2}_?\d+\.\d+)(?:_\d+)?$", tip_name)
        if m_new:
            label = f"{m_new.group(1).replace('_', ' ')} ({m_new.group(2)})"
        else:
            label = tip_name.replace("_", " ")
        return False, "#334155", 3.0, label

    row_height = 22
    margin_top = 40
    margin_bottom = 60
    margin_left = 40
    margin_right = 340
    tree_width = 680
    svg_width = margin_left + tree_width + margin_right
    svg_height = margin_top + n_leaves * row_height + margin_bottom

    def to_px_x(val):
        return margin_left + (val / max_x) * tree_width

    def to_px_y(val):
        return margin_top + val * row_height

    svg_id = f"mitoconda_phylo_svg_{phylo_target}_{pident}_{gene}"
    svg_lines = [
        f'<svg id="{svg_id}" data-base-w="{svg_width}" data-base-h="{svg_height}" data-zoom="1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_width} {svg_height}" width="{svg_width}" height="{svg_height}" style="background:#ffffff; font-family: system-ui, -apple-system, sans-serif; display:block;">',
        f'<rect width="{svg_width}" height="{svg_height}" fill="#ffffff" />',
    ]

    for parent in pruned_tree.find_clades(order="level"):
        px = to_px_x(x_coords[id(parent)])
        if parent.clades:
            child_ys = [to_px_y(y_coords[id(c)]) for c in parent.clades]
            min_cy, max_cy = min(child_ys), max(child_ys)
            svg_lines.append(
                f'<line x1="{px:.2f}" y1="{min_cy:.2f}" x2="{px:.2f}" y2="{max_cy:.2f}" stroke="#64748b" stroke-width="1.4" />'
            )
            for child in parent.clades:
                cx = to_px_x(x_coords[id(child)])
                cy = to_px_y(y_coords[id(child)])
                svg_lines.append(
                    f'<line x1="{px:.2f}" y1="{cy:.2f}" x2="{cx:.2f}" y2="{cy:.2f}" stroke="#64748b" stroke-width="1.4" />'
                )
                if child.clades and child.confidence is not None:
                    try:
                        conf = float(child.confidence)
                        conf_pct = int(round(conf * 100)) if conf <= 1.0 else int(round(conf))
                        if conf_pct >= 70:
                            mid_x = (px + cx) / 2.0
                            svg_lines.append(
                                f'<text x="{mid_x:.2f}" y="{cy - 4:.2f}" font-size="9" fill="#64748b" text-anchor="middle">{conf_pct}</text>'
                            )
                    except ValueError:
                        pass

    for leaf in leaves:
        lx = to_px_x(x_coords[id(leaf)])
        ly = to_px_y(y_coords[id(leaf)])
        is_sample, col, rad, label_str = get_tip_style(leaf.name or "")
        esc_label = html.escape(label_str)
        if is_sample:
            svg_lines.append(
                f'<circle cx="{lx:.2f}" cy="{ly:.2f}" r="{rad:.1f}" fill="{col}" stroke="#ffffff" stroke-width="1.2" />'
            )
            svg_lines.append(
                f'<text x="{lx + 10:.2f}" y="{ly + 4:.2f}" font-size="12" font-weight="700" fill="{col}">{esc_label}</text>'
            )
        else:
            svg_lines.append(
                f'<circle cx="{lx:.2f}" cy="{ly:.2f}" r="{rad:.1f}" fill="#94a3b8" />'
            )
            svg_lines.append(
                f'<text x="{lx + 10:.2f}" y="{ly + 4:.2f}" font-size="11.5" font-style="italic" fill="{col}">{esc_label}</text>'
            )

    raw_scale = max_x / 5.0
    magnitude = 10 ** math.floor(math.log10(raw_scale)) if raw_scale > 0 else 0.01
    scale_val = round(raw_scale / magnitude) * magnitude
    if scale_val <= 0:
        scale_val = 0.05
    scale_px = (scale_val / max_x) * tree_width
    sb_x1 = margin_left
    sb_x2 = margin_left + scale_px
    sb_y = svg_height - 25
    svg_lines.extend([
        f'<line x1="{sb_x1:.2f}" y1="{sb_y:.2f}" x2="{sb_x2:.2f}" y2="{sb_y:.2f}" stroke="#1e293b" stroke-width="2" />',
        f'<line x1="{sb_x1:.2f}" y1="{sb_y - 4:.2f}" x2="{sb_x1:.2f}" y2="{sb_y + 4:.2f}" stroke="#1e293b" stroke-width="2" />',
        f'<line x1="{sb_x2:.2f}" y1="{sb_y - 4:.2f}" x2="{sb_x2:.2f}" y2="{sb_y + 4:.2f}" stroke="#1e293b" stroke-width="2" />',
        f'<text x="{(sb_x1 + sb_x2) / 2:.2f}" y="{sb_y + 16:.2f}" font-size="11" fill="#1e293b" text-anchor="middle">{scale_val:g} substitutions/site</text>',
        "</svg>",
    ])

    with open(svg_file, "w") as sf:
        sf.write("\n".join(svg_lines))

    branch_x = []
    branch_y = []
    support_x = []
    support_y = []
    support_text = []
    support_hover = []

    for parent in pruned_tree.find_clades(order="level"):
        px = float(x_coords[id(parent)])
        if parent.clades:
            child_ys = [float(y_coords[id(c)]) for c in parent.clades]
            min_cy, max_cy = min(child_ys), max(child_ys)
            branch_x.extend([px, px, None])
            branch_y.extend([min_cy, max_cy, None])
            for child in parent.clades:
                cx = float(x_coords[id(child)])
                cy = float(y_coords[id(child)])
                branch_x.extend([px, cx, None])
                branch_y.extend([cy, cy, None])
                if child.clades and child.confidence is not None:
                    try:
                        conf = float(child.confidence)
                        conf_pct = int(round(conf * 100)) if conf <= 1.0 else int(round(conf))
                        if conf_pct >= 70:
                            support_x.append((px + cx) / 2.0)
                            support_y.append(cy)
                            support_text.append(f"{conf_pct}%")
                            support_hover.append(f"FastTree Support: {conf_pct}%<br>Branch Length: {cx - px:.5f}")
                    except ValueError:
                        pass

    plotly_traces = [
        {
            "type": "scatter",
            "mode": "lines",
            "name": "Tree Branches",
            "x": branch_x,
            "y": branch_y,
            "line": {"color": "#64748b", "width": 1.5},
            "hoverinfo": "skip",
            "showlegend": False,
        }
    ]

    if support_x:
        plotly_traces.append({
            "type": "scatter",
            "mode": "markers+text",
            "name": "Support >= 70%",
            "x": support_x,
            "y": support_y,
            "text": support_text,
            "textposition": "top center",
            "textfont": {"size": 9, "color": "#475569"},
            "marker": {"color": "#475569", "size": 4, "symbol": "diamond"},
            "hovertext": support_hover,
            "hoverinfo": "text",
        })

    ref_x, ref_y, ref_text, ref_hover = [], [], [], []
    bc_traces_data = defaultdict(lambda: {"x": [], "y": [], "text": [], "hover": [], "sizes": [], "color": "#2563eb"})

    for leaf in leaves:
        lx = float(x_coords[id(leaf)])
        ly = float(y_coords[id(leaf)])
        tname = leaf.name or ""
        is_sample, col, rad, label_str = get_tip_style(tname)
        own_bl = float(leaf.branch_length) if leaf.branch_length is not None else 0.0
        if is_sample:
            matched_sample = "Sample"
            for s in sorted(samples):
                short_s = re.sub(r"_junto$", "", s)
                if tname.startswith(f"{s}_") or tname.startswith(f"{short_s}_"):
                    matched_sample = short_s
                    break
            m_sz = re.search(r"size[_=]?(\d+)", tname)
            reads = int(m_sz.group(1)) if m_sz else 1
            near_ref, near_dist = nearest_ref_for_cluster.get(id(leaf), ("N/A", 0.0))
            grp = bc_traces_data[matched_sample]
            grp["x"].append(lx)
            grp["y"].append(ly)
            grp["text"].append(f"  <b>{label_str}</b>")
            grp["sizes"].append(max(8.0, min(16.0, rad * 2.1)))
            grp["color"] = col
            grp["hover"].append(
                f"<b>{label_str}</b><br>"
                f"Sample: {matched_sample}<br>"
                f"Reads: {reads}<br>"
                f"Closest Reference: <i>{near_ref}</i><br>"
                f"Patristic Distance: {near_dist:.5f} subs/site<br>"
                f"Terminal Branch: {own_bl:.5f} subs/site"
            )
        else:
            ref_x.append(lx)
            ref_y.append(ly)
            ref_text.append(f"  <i>{label_str}</i>")
            ref_hover.append(
                f"<b>Reference Taxon</b><br>"
                f"Species: <i>{label_str}</i><br>"
                f"Root Distance: {lx:.5f} subs/site<br>"
                f"Terminal Branch: {own_bl:.5f} subs/site"
            )

    if ref_x:
        plotly_traces.append({
            "type": "scatter",
            "mode": "markers+text",
            "name": "Reference Taxa (NCBI Mito)",
            "x": ref_x,
            "y": ref_y,
            "text": ref_text,
            "textposition": "middle right",
            "textfont": {"size": 11, "color": "#334155"},
            "marker": {"color": "#94a3b8", "size": 6},
            "hovertext": ref_hover,
            "hoverinfo": "text",
            "cliponaxis": False,
        })

    for s_name in sorted(bc_traces_data.keys()):
        grp = bc_traces_data[s_name]
        plotly_traces.append({
            "type": "scatter",
            "mode": "markers+text",
            "name": s_name,
            "x": grp["x"],
            "y": grp["y"],
            "text": grp["text"],
            "textposition": "middle right",
            "textfont": {"size": 12, "color": grp["color"]},
            "marker": {
                "color": grp["color"],
                "size": grp["sizes"],
                "line": {"color": "#ffffff", "width": 1.2},
            },
            "hovertext": grp["hover"],
            "hoverinfo": "text",
            "cliponaxis": False,
        })

    plot_height = max(650, min(1400, n_leaves * 18 + 120))
    plotly_layout = {
        "title": {
            "text": f"Integrated Phylogenetic Tree ({len(cluster_tips)} Sample Clusters & {len(leaves) - len(cluster_tips)} Closest Reference Taxa <= 5 Nodes)",
            "font": {"size": 15},
        },
        "height": plot_height,
        "margin": {"l": 30, "r": 260, "t": 50, "b": 55},
        "hovermode": "closest",
        "dragmode": "zoom",
        "plot_bgcolor": "#ffffff",
        "paper_bgcolor": "#ffffff",
        "xaxis": {
            "title": {"text": "Genetic Distance (substitutions / site)"},
            "range": [-max_x * 0.02, max_x * 1.48],
            "showgrid": True,
            "gridcolor": "#f1f5f9",
            "zeroline": False,
        },
        "yaxis": {
            "autorange": "reversed",
            "showticklabels": False,
            "showgrid": False,
            "zeroline": False,
        },
        "legend": {
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.01,
            "xanchor": "left",
            "x": 0,
        },
    }
    plotly_config = {
        "responsive": True,
        "displaylogo": False,
        "scrollZoom": True,
        "toImageButtonOptions": {
            "format": "svg",
            "filename": f"All_Barcodes_{phylo_target}_{pident}_{gene}_tree",
        },
    }

    plotly_div_id = f"mqc_plotly_phylo_{phylo_target}_{pident}_{gene}"
    traces_json = json.dumps(plotly_traces)
    layout_json = json.dumps(plotly_layout)
    config_json = json.dumps(plotly_config)

    html_report = f"""<!--
id: 'phylo_tree_visualization'
section_name: '3.1 Integrated Phylogenetic Tree (All Barcodes vs. Closest References)'
description: '<b>Methodology:</b> Consensus clusters across all samples matching the target clade are trimmed to the aligned 12S region, oriented to the forward strand, and inserted into the taxonomically balanced NCBI mitochondrial reference alignment via <code>MAFFT --addfragments --keeplength</code>. Flanking reference-only regions are removed, gappy columns are trimmed with <code>trimAl (-gt 0.7)</code>, and a Maximum-Likelihood tree is inferred with <code>FastTree (-nt -gtr -gamma)</code>. The full tree is midpoint-rooted prior to pruning to display all sample clusters alongside reference branches within <b>5 internal nodes of topological distance</b> (preserving exact root-to-tip branch lengths). <b>Interactive Controls:</b> Colored circles represent sample clusters (circle size is proportional to log<sub>10</sub> read abundance); gray nodes represent NCBI reference mitogenomes; diamond markers show FastTree local support values &ge;70%. Drag a box to zoom into any clade, use the scroll wheel to zoom, double-click to reset, hover over tips for exact patristic distances, or click samples in the top legend to toggle visibility.'
plot_type: 'html'
-->
<div class="mqc-custom-content-plot" style="border:1px solid #e2e8f0; border-radius:8px; padding:10px; background:#ffffff;">
    <div id="{plotly_div_id}" style="width:100%; height:{plot_height}px;"></div>
</div>
<script type="text/javascript">
(function renderMitoCondaPhyloTree() {{
    if (typeof Plotly === "undefined") {{
        if (!window._mitoconda_loading_plotly) {{
            window._mitoconda_loading_plotly = true;
            var s = document.createElement("script");
            s.src = "https://cdn.plot.ly/plotly-2.35.2.min.js";
            document.head.appendChild(s);
        }}
        setTimeout(renderMitoCondaPhyloTree, 100);
        return;
    }}
    var target = document.getElementById("{plotly_div_id}");
    if (!target) return;
    var traces = {traces_json};
    var layout = {layout_json};
    var config = {config_json};
    Plotly.newPlot(target, traces, layout, config);
}})();
</script>
"""
    with open(mqc_html, "w") as hf:
        hf.write(html_report)

