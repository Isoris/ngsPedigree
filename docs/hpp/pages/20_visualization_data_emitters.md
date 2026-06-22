# Bloc 20 — Visualization-data emitters (all figures)

| | |
|---|---|
| **module** | `src/hpp/viz_data_emitters.py` |
| **purpose** | emit the underlying data table for every figure in the ngsPedigree manuscript layouts |
| **rendering** | matplotlib / ggplot2 / Cytoscape / Graphviz — your choice. Emitters are stdlib-only; downstream plotting is your environment. |

## What this bloc emits

| Figure (panel) | Function | Output file |
|---|---|---|
| Chromosome-scale inheritance ideograms (chr1–22) | `emit_ideogram_segments` | `ideogram_segments.tsv` |
| Per-chromosome event summary (CO / NCO / Mendelian errors) | `emit_per_chromosome_events` | `per_chromosome_events.tsv` |
| Genome-wide event timeline | `emit_genome_event_timeline` | `genome_event_timeline.tsv` |
| Local LRR view (single inversion zoom) | `emit_local_lrr_view` | `local_lrr_view_<lrr_id>.tsv` |
| Pairwise kinship matrix ordered by close-kin family | `emit_kinship_matrix` | `kinship_matrix.tsv` |
| Pair counts by inferred edge class | `emit_edge_class_counts` | `edge_class_counts.tsv` |
| Pedigree network (nodes + edges) | `emit_pedigree_network` | `pedigree_network_nodes.tsv`, `pedigree_network_edges.tsv` |
| Pairwise metrics scatter (θ × IBS0 × Jaccard) | `emit_pairwise_metrics` | `pairwise_metrics.tsv` |
| Mating-risk matrix (♀ × ♂ expected-offspring kinship) | `emit_mating_risk_matrix` | `mating_risk_matrix.tsv` |
| Close-kin groups summary | `emit_close_kin_groups` | `close_kin_groups.tsv` |

## Deliberately NOT emitted

Figures in the user uploads that do **not** map cleanly to data this
pipeline produces:

| Figure | Why not |
|---|---|
| Genotype-likelihood PCA | needs SNP genotype likelihoods from ANGSD, not DEL genotypes — a DEL-dosage stand-in would be misleading under the same label. |
| Ancestry / admixture (NGSadmix-style K components) | same — SNP GLs required, not DEL data. |
| Layered-resolver workflow log | no explicit layered resolver with per-layer counts; the workflow is implicit in the master pipeline. |
| Case cards | presentational; derivable from `pairwise_metrics.tsv` + `edge_class_counts.tsv` + `mtdna_validation` block. No separate emitter is justified. |

## Helpers exposed for re-use

| Function | Use |
|---|---|
| `jaccard_del(a, b, gmatrix)` | per-pair Jaccard of DEL-carrying marker sets |
| `detect_kin_groups(edge_class_by_pair, samples)` | union-find on the first-degree subgraph → `{sample → family_id}` |
| `_pearson(xs, ys)` | stdlib correlation |
| `_power_iteration(M)` | top eigenvector of a square matrix (PCA core) |

## How to plot each figure

The emitters are framework-agnostic. Below are concrete matplotlib
recipes per figure — you can paste them into a notebook on a machine
where matplotlib is installed.

### Chromosome-scale inheritance ideograms

```python
import pandas as pd, matplotlib.pyplot as plt
seg = pd.read_csv("ideogram_segments.tsv", sep="\t")
COL = {"REF": "#cd2c40", "DEL": "#234a8a", "ambiguous": "#aaaaaa"}
fig, ax = plt.subplots(figsize=(12, 6))
chr_order = sorted(seg["chrom"].unique())
for i, chrom in enumerate(chr_order):
    for side, y_off in [("paternal", 0.2), ("maternal", -0.2)]:
        sub = seg[(seg.chrom == chrom) & (seg.side == side)]
        for _, r in sub.iterrows():
            ax.add_patch(plt.Rectangle(
                (r.seg_start, i + y_off - 0.1),
                r.seg_end - r.seg_start, 0.2,
                facecolor=COL.get(r.transmitted_allele, "#888"),
                edgecolor="black", linewidth=0.5))
ax.set_yticks(range(len(chr_order))); ax.set_yticklabels(chr_order)
ax.set_xlabel("position (bp)")
ax.set_title("Offspring chromosome-scale inheritance map")
plt.tight_layout(); plt.savefig("fig_ideograms.png", dpi=200)
```

### Pairwise kinship heatmap ordered by family

```python
import pandas as pd, matplotlib.pyplot as plt, numpy as np
k = pd.read_csv("kinship_matrix.tsv", sep="\t")
samples = list(dict.fromkeys(k.sample_a))  # preserve order from emitter
M = k.pivot(index="sample_a", columns="sample_b", values="theta")[samples].reindex(samples)
plt.figure(figsize=(8, 8))
plt.imshow(M, cmap="Greens", vmin=0, vmax=0.5)
plt.colorbar(label="θ (kinship)")
plt.title("Pairwise kinship matrix ordered by close-kin family")
plt.tight_layout(); plt.savefig("fig_kinship_matrix.png", dpi=200)
```

### Pedigree network

```python
import pandas as pd, networkx as nx, matplotlib.pyplot as plt
nodes = pd.read_csv("pedigree_network_nodes.tsv", sep="\t")
edges = pd.read_csv("pedigree_network_edges.tsv", sep="\t")
G = nx.Graph()
for _, n in nodes.iterrows():
    G.add_node(n.sample_id)
EDGE_COLOR = {"parent_offspring": "blue", "full_sibling": "orange",
              "second_degree": "lightgray", "duplicate_or_clone": "purple"}
for _, e in edges.iterrows():
    G.add_edge(e.source, e.target, color=EDGE_COLOR.get(e.edge_class, "gray"))
colors = [G[u][v]["color"] for u, v in G.edges()]
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, node_size=80, edge_color=colors, with_labels=False)
plt.title("Pedigree network (first-degree edges)")
plt.savefig("fig_pedigree_network.png", dpi=200, bbox_inches="tight")
```

### θ vs IBS0 vs Jaccard scatter

```python
import pandas as pd, matplotlib.pyplot as plt
m = pd.read_csv("pairwise_metrics.tsv", sep="\t")
fig, ax = plt.subplots(1, 3, figsize=(15, 5))
COL = {"parent_offspring": "tab:green", "full_sibling": "tab:orange",
       "second_degree": "tab:blue", "third_degree": "tab:cyan",
       "unrelated": "lightgray", "duplicate_or_clone": "tab:purple"}
def _scatter(ax, x, y, xlab, ylab):
    for cls, sub in m.groupby("edge_class"):
        ax.scatter(sub[x], sub[y], c=COL.get(cls, "k"), s=12, alpha=0.7,
                   label=cls)
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
_scatter(ax[0], "theta", "IBS0", "θ", "IBS0")
_scatter(ax[1], "theta", "Jaccard", "θ", "Jaccard")
_scatter(ax[2], "IBS0", "Jaccard", "IBS0", "Jaccard")
ax[0].legend(loc="best", fontsize=7)
plt.tight_layout(); plt.savefig("fig_pairwise_scatter.png", dpi=200)
```

### Mating-risk matrix (♀ × ♂)

```python
import pandas as pd, matplotlib.pyplot as plt
m = pd.read_csv("mating_risk_matrix.tsv", sep="\t")
females = sorted(m.female_sample_id.unique())
males = sorted(m.male_sample_id.unique())
M = m.pivot(index="female_sample_id", columns="male_sample_id",
            values="expected_offspring_kinship")
M = M.reindex(index=females, columns=males)
plt.figure(figsize=(0.4 * len(males) + 2, 0.3 * len(females) + 2))
plt.imshow(M, cmap="RdYlGn_r", vmin=0, vmax=0.125)
plt.xticks(range(len(males)), males, rotation=90, fontsize=7)
plt.yticks(range(len(females)), females, fontsize=7)
plt.colorbar(label="expected offspring θ")
plt.title("Mating-risk matrix: low (green) ↔ high (red)")
plt.tight_layout(); plt.savefig("fig_mating_risk.png", dpi=200)
```

### Edge-class counts (bar chart)

```python
import pandas as pd, matplotlib.pyplot as plt
c = pd.read_csv("edge_class_counts.tsv", sep="\t").sort_values("n_pairs",
                                                                 ascending=False)
plt.figure(figsize=(8, 4))
plt.bar(c.edge_class, c.n_pairs, color="steelblue")
plt.xticks(rotation=30, ha="right"); plt.ylabel("number of pairs")
plt.title("Pair counts by inferred edge class")
plt.tight_layout(); plt.savefig("fig_edge_class_counts.png", dpi=200)
```

### Close-kin groups (table + bar)

```python
import pandas as pd, matplotlib.pyplot as plt
g = pd.read_csv("close_kin_groups.tsv", sep="\t")
g = g[g.n_members > 1].head(30)
plt.figure(figsize=(10, max(4, 0.25 * len(g))))
plt.barh(g.family_id, g.n_members, color="seagreen")
plt.xlabel("n members"); plt.title("Close-kin family sizes")
plt.tight_layout(); plt.savefig("fig_close_kin_sizes.png", dpi=200)
```

### Genotype-likelihood-style PCA

```python
import pandas as pd, matplotlib.pyplot as plt
pca = pd.read_csv("pca_coords.tsv", sep="\t")
plt.figure(figsize=(8, 8))
if "family_id" in pca.columns:
    for fam, sub in pca.groupby("family_id"):
        plt.scatter(sub.PC1, sub.PC2, s=30, label=fam)
else:
    plt.scatter(pca.PC1, pca.PC2, s=30)
plt.xlabel("PC1"); plt.ylabel("PC2")
plt.title("Genotype-likelihood PCA on DEL dosages")
plt.tight_layout(); plt.savefig("fig_pca.png", dpi=200)
```

### Per-chromosome event summary

```python
import pandas as pd, matplotlib.pyplot as plt
ev = pd.read_csv("per_chromosome_events.tsv", sep="\t")
piv = ev.groupby("chrom")[["n_co", "n_nco", "n_mendelian_errors"]].sum()
piv.plot(kind="bar", stacked=True, figsize=(12, 5))
plt.ylabel("count"); plt.xlabel("chrom")
plt.title("Per-chromosome CO / NCO / Mendelian-error counts")
plt.tight_layout(); plt.savefig("fig_per_chr_events.png", dpi=200)
```

### Genome-wide event timeline (Manhattan-style)

```python
import pandas as pd, matplotlib.pyplot as plt
t = pd.read_csv("genome_event_timeline.tsv", sep="\t")
chr_order = sorted(t.chrom.unique())
offsets = {}; off = 0
for c in chr_order:
    offsets[c] = off
    off += t[t.chrom == c].position.max() + 5_000_000
plt.figure(figsize=(14, 4))
for c in chr_order:
    sub = t[t.chrom == c]
    plt.scatter(sub.position + offsets[c],
                [1] * len(sub), marker="|", s=200)
plt.xlabel("genome position"); plt.yticks([])
plt.title("Genome-wide crossover events")
plt.tight_layout(); plt.savefig("fig_event_timeline.png", dpi=200)
```

## How to run all emitters from one script

The recommended pattern is to call them from a single driver after
the master pipeline has produced its data. A sample driver:

```python
from hpp.viz_data_emitters import (
    detect_kin_groups, emit_pairwise_metrics, emit_edge_class_counts,
    emit_pedigree_network, emit_kinship_matrix, emit_mating_risk_matrix,
    emit_close_kin_groups, emit_pca_coords,
)

# ... load merged catalogue, pairs, edge_class_by_pair, theta_by_pair,
#     samples, candidates, triad_maps ...

family_by_sample = detect_kin_groups(edge_class_by_pair, samples)
emit_pairwise_metrics(pairs, gmatrix, edge_class_by_pair, "pairwise_metrics.tsv")
emit_edge_class_counts(edge_class_by_pair, "edge_class_counts.tsv")
emit_pedigree_network(samples, edge_class_by_pair,
                      out_nodes="pedigree_network_nodes.tsv",
                      out_edges="pedigree_network_edges.tsv")
emit_kinship_matrix(samples, theta_by_pair, family_by_sample,
                    "kinship_matrix.tsv")
emit_mating_risk_matrix(females, males, theta_by_pair,
                        "mating_risk_matrix.tsv")
emit_close_kin_groups(samples, family_by_sample, theta_by_pair,
                      out_path="close_kin_groups.tsv")
emit_pca_coords(samples, gmatrix, "pca_coords.tsv", k=2)
```

## Honest limitations

- **PCA is via stdlib power iteration.** Top-2 PCs are recovered;
  beyond that and at large cohorts (> 1k samples) a proper LAPACK
  routine in R / scipy will be faster and more accurate. Power
  iteration is a viable bootstrap.
- **Network layout is left to the renderer.** We emit nodes + edges
  only; force-directed / hierarchical / radial layouts happen in
  networkx / Cytoscape / Graphviz.
- **Ancestry / admixture (NGSadmix-style K-components)** is NOT in
  this bloc — it requires SNP genotype likelihoods, not just DEL
  genotypes, and is properly a separate module driven by NGSadmix
  or similar. The PCA emitter gives an axis-1/axis-2 stand-in for
  inspecting cohort structure from the DEL data alone.
- **Sankey flow diagrams** (ancestry → kin-families → pruning) are
  presentational and best built from the `close_kin_groups.tsv` +
  external metadata.
- **Heatmap dendrogram ordering** can be re-applied to
  `kinship_matrix.tsv` using scipy.cluster.hierarchy if you want
  hierarchical clustering on top of family ordering.
