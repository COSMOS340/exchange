"""Full-cell plots of saved CRC models. This script never fits or samples cells."""
import argparse
import json
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/02_表格/CRC"
OUT = ROOT / "results/01_图形/main"
LOG = ROOT / "results/04_日志和检查"
BANDS = ["Low", "Intermediate", "High"]
COLORS = {"Low": "#3182BD", "Intermediate": "#BDBDBD", "High": "#D94B4B"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.facecolor": "white"})


def aligned(filename, metadata):
    data = pd.read_csv(TAB / filename)
    assert data.cell_id.is_unique and set(data.cell_id) == set(metadata.cell_id)
    return data.set_index("cell_id").loc[metadata.cell_id].reset_index()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slingshot", action="store_true")
    parser.add_argument("--only", nargs="*", help="Re-render only named figures after a visual repair")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(ROOT / "work/CRC/metadata.csv.gz", low_memory=False)
    assert len(metadata) == 150917 and metadata.cell_id.is_unique
    main_cells = aligned("Monocle2_primary_cells.csv.gz", metadata)
    alt_cells = aligned("Monocle2_alternative_cells.csv.gz", metadata)
    xy = main_cells[["DDRTree_1", "DDRTree_2"]].to_numpy()
    assert np.isfinite(xy).all()
    assert np.array_equal(xy, alt_cells[["DDRTree_1", "DDRTree_2"]].to_numpy())
    centers = pd.read_csv(TAB / "Monocle2_centers.csv").set_index("center")
    centers.index = centers.index.astype(str)
    edges = pd.read_csv(TAB / "Monocle2_tree_edges.csv")
    roots = json.loads((TAB / "Monocle2_roots.json").read_text())
    permutation = np.random.default_rng(20261020).permutation(len(metadata))
    assert np.array_equal(np.sort(permutation), np.arange(len(metadata)))
    span = np.ptp(xy, axis=0)
    limits = [(float(xy[:, j].min() - .04 * span[j]),
               float(xy[:, j].max() + .04 * span[j])) for j in range(2)]
    manifest = []
    band_counts = metadata.CNV_relative_band.value_counts().reindex(BANDS)
    source_cells = "02_表格/CRC/Monocle2_primary_cells.csv.gz"

    def save(fig, name, question, finding, sources, unit="cell-level display"):
        png, pdf = OUT / (name + ".png"), OUT / (name + ".pdf")
        if not args.only or name in args.only:
            # Solve the layout once before switching PNG/PDF backends. A
            # second constrained-layout pass can otherwise clip table artists.
            fig.canvas.draw()
            fig.set_layout_engine("none")
            fig.savefig(png, dpi=180, bbox_inches="tight")
            fig.savefig(pdf, dpi=180, bbox_inches="tight")
        else:
            assert png.exists() and pdf.exists(), f"Missing skipped figure: {name}"
        plt.close(fig)
        manifest.append(dict(name=name, png=str(png.relative_to(ROOT / "results")),
            pdf=str(pdf.relative_to(ROOT / "results")), question=question,
            finding=finding, sources=sources, statistical_unit=unit,
            fitted_cells=len(metadata), randomly_omitted_cells=0,
            visual_review="pending"))

    def tree(ax):
        for a, b in edges[["from", "to"]].itertuples(index=False, name=None):
            segment = centers.loc[[str(a), str(b)], ["DDRTree_1", "DDRTree_2"]].to_numpy()
            ax.plot(*segment.T, color="#252525", lw=.65, alpha=.65, zorder=3)

    def clean(ax, title):
        ax.set(title=title, xlabel="DDRTree 1", ylabel="DDRTree 2",
               xlim=limits[0], ylim=limits[1])
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(labelsize=7)

    def bands(ax, selected=None, title="Relative RNA CNV", legend=True):
        keep = np.ones(len(metadata), dtype=bool) if selected is None else np.asarray(selected)
        if selected is not None:
            ax.scatter(*xy.T, s=.3, c="#E8E8E8", marker="o", linewidths=0, rasterized=True)
        order = permutation[keep[permutation]]
        ax.scatter(*xy[order].T, s=.7, c=metadata.CNV_relative_band.map(COLORS).to_numpy()[order],
                   marker="o", linewidths=0, alpha=.75, rasterized=True)
        tree(ax); clean(ax, title)
        if legend:
            handles = [Line2D([], [], marker="o", linestyle="", color=COLORS[b],
                       label=f"{b}: n={band_counts[b]:,}", markersize=5) for b in BANDS]
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1),
                      frameon=False, fontsize=8)

    def numeric(fig, ax, values, title, label, cmap="viridis", domain=None):
        v = np.asarray(values, float); assert np.isfinite(v).all()
        lo, hi = (float(v.min()), float(v.max())) if domain is None else domain
        artist = ax.scatter(*xy[permutation].T, c=v[permutation], s=.7, cmap=cmap,
            vmin=lo, vmax=hi if hi > lo else lo + 1, marker="o", linewidths=0, rasterized=True)
        tree(ax); clean(ax, title)
        fig.colorbar(artist, ax=ax, shrink=.65, label=label, pad=.02)

    def category(ax, values, title):
        v = pd.Series(values).fillna("Missing").astype(str).to_numpy()
        levels = sorted(set(v), key=lambda x: (len(x), x))
        cmap = plt.get_cmap("tab20" if len(levels) <= 20 else "turbo")
        palette = {level: cmap(i / max(1, len(levels)-1)) for i, level in enumerate(levels)}
        ax.scatter(*xy[permutation].T, c=[palette[x] for x in v[permutation]], s=.65,
                   marker="o", linewidths=0, rasterized=True)
        tree(ax); clean(ax, title)
        if len(levels) <= 30:
            handles = [Line2D([], [], marker="o", linestyle="", markersize=4,
                       color=palette[x], label=textwrap.fill(f"{x}: {(v == x).sum():,}", 24))
                       for x in levels]
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
                      frameon=False, fontsize=6)
        else:
            mapping = pd.DataFrame({"label": levels, "color_RGBA": [str(palette[x]) for x in levels]})
            mapping.to_csv(TAB / (title.replace(" ", "_") + "_color_key.csv"), index=False)
            ax.text(0, -.16, "Full color key is in the source table", transform=ax.transAxes, fontsize=7)

    description = metadata[["cell_id", "dataset", "sample_id", "specimen_id", "CNV_relative_band"]].copy()
    description["pseudotime"] = main_cells.PT_1
    pt_stats = description.groupby("CNV_relative_band").pseudotime.agg(
        cells="size", median="median", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75)).reindex(BANDS)
    pt_stats.to_csv(TAB / "CNV_pseudotime_description.csv")
    finding = "；".join(f"{b}组{int(pt_stats.loc[b, 'cells']):,}细胞，拟时序中位数{pt_stats.loc[b, 'median']:.3g}" for b in BANDS)
    fig, ax = plt.subplots(figsize=(9, 7), layout="constrained")
    bands(ax, title="CRC-origin epithelial cells | relative RNA CNV")
    save(fig, "01_Monocle2_CNV_groups", "既有CNV三组在同一树上的分布如何？", finding,
         [source_cells, "02_表格/CRC/CNV_pseudotime_description.csv"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), layout="constrained")
    pt_domain = (0., float(max(main_cells.PT_1.max(), alt_cells.PT_1.max())))
    numeric(fig, axes[0], main_cells.PT_1, "Computational-reference pseudotime", "Weighted tree distance (a.u.)", domain=pt_domain)
    numeric(fig, axes[1], metadata.cnv_abs_log2, "Continuous relative RNA CNV", "cnv_abs_log2 (mean absolute log2)")
    save(fig, "02_Monocle2_PT_continuous_CNV", "拟时序和连续CNV在共同表达轨迹上的位置如何？",
         finding + "。共同坐标不证明CNV随真实时间变化。", [source_cells])

    fig, axes = plt.subplots(1, 3, figsize=(16, 6), layout="constrained")
    for ax, band in zip(axes, BANDS):
        bands(ax, metadata.CNV_relative_band.eq(band), f"{band} | n={band_counts[band]:,}", legend=False)
    save(fig, "03_Monocle2_CNV_facets", "三组在完全相同坐标和完整树下有何覆盖？",
         "三面板分别显示既有Low/Intermediate/High组，其他细胞为灰色背景；各面板保留全部拟合细胞的背景和同一棵树。",
         [source_cells])

    for number, field, title in [(4, "State", "Monocle2 State"), (5, "coarse_state", "Inherited epithelial cluster"),
                                  (6, "broad_identity", "Inherited broad identity"), (7, "dataset", "Source cohort"),
                                  (8, "Phase", "Inherited cell-cycle phase")]:
        values = main_cells.State if field == "State" else metadata[field]
        fig, ax = plt.subplots(figsize=(10, 8), layout="constrained")
        category(ax, values, title)
        counts = pd.Series(values).value_counts().rename_axis(field).reset_index(name="cells")
        counts.to_csv(TAB / f"{field}_coverage.csv", index=False)
        save(fig, f"{number:02d}_Monocle2_{field}", f"同一树按{title}着色有何分布？",
             f"共{len(counts)}个{title}类别；最多类别为{counts.iloc[0,0]}，含{counts.iloc[0,1]:,}细胞。State是本次树区段，继承簇和宽身份是上游标签，不能互换。",
             [source_cells, f"02_表格/CRC/{field}_coverage.csv"])

    def composition(frame, column, prefix, title):
        counts = pd.crosstab(frame[column], frame.CNV_relative_band).reindex(columns=BANDS, fill_value=0)
        pooled = counts.div(counts.sum(axis=1), axis=0)
        keys = [column, "dataset", "specimen_id"]
        units = frame.groupby(keys + ["CNV_relative_band"], observed=True).size().unstack(fill_value=0).reindex(columns=BANDS, fill_value=0)
        equal = units.div(units.sum(axis=1), axis=0).groupby(level=0).mean().reindex(pooled.index)
        coverage = frame.groupby(column).agg(cells=("cell_id", "size"),
            specimens=("specimen_id", "nunique"), libraries=("sample_id", "nunique"),
            reported_labels=("subject_key", "nunique")).reindex(pooled.index)
        counts.to_csv(TAB / f"{prefix}_CNV_counts.csv")
        pooled.to_csv(TAB / f"{prefix}_CNV_pooled_fraction.csv")
        equal.to_csv(TAB / f"{prefix}_CNV_specimen_equal_fraction.csv")
        coverage.to_csv(TAB / f"{prefix}_sample_coverage.csv")
        units.div(units.sum(axis=1), axis=0).to_csv(TAB / f"{prefix}_CNV_per_specimen_fraction.csv")
        fig, axes = plt.subplots(1, 3, figsize=(13, max(5, .26*len(pooled))),
                                 layout="constrained", gridspec_kw={"width_ratios": [1, 1, .85]})
        for ax, values, heading in zip(axes[:2], [pooled, equal], ["Pooled-cell proportion", "Equal-specimen mean proportion"]):
            image = ax.imshow(values.to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="Blues")
            ax.set(xticks=range(3), xticklabels=BANDS, yticks=range(len(values)),
                   yticklabels=values.index.astype(str), title=heading)
            ax.tick_params(labelsize=7)
            fig.colorbar(image, ax=ax, shrink=.6, label="Fraction (0-1)")
        axes[2].axis("off")
        axes[2].set_title("Contributing coverage")
        coverage_table = axes[2].table(cellText=[[f"{r.cells:,}", f"{r.specimens:,}"] for r in coverage.itertuples()],
                      rowLabels=coverage.index.astype(str), colLabels=["Cells", "Specimens"], loc="center", cellLoc="center")
        coverage_table.auto_set_font_size(False)
        coverage_table.set_fontsize(9)
        coverage_table.scale(1, 1.35)
        fig.suptitle(title)
        save(fig, f"09_{prefix}_CNV_composition", f"{title}的CNV组成和标本覆盖如何？",
             f"展示{len(pooled)}个区段/路径的细胞合并比例与标本等权比例；样本覆盖见右侧及完整源表。共享树干会重复出现在不同路径，路径计数不能累加为独立细胞。",
             [f"02_表格/CRC/{prefix}_CNV_pooled_fraction.csv",
              f"02_表格/CRC/{prefix}_CNV_specimen_equal_fraction.csv", f"02_表格/CRC/{prefix}_sample_coverage.csv"],
             unit="pooled cells and equal contributing specimens, shown separately")

    state_frame = metadata.copy(); state_frame["State"] = main_cells.State.astype(str)
    composition(state_frame, "State", "Monocle2_State", "CNV composition by Monocle2 State")
    for orientation in ["primary", "alternative"]:
        path_cells = aligned(f"Monocle2_{orientation}_path_cells.csv.gz", metadata)
        parts = []
        for field in [c for c in path_cells if c.startswith("PT_")]:
            number = field.split("_")[1]
            keep = path_cells[field].notna() & path_cells["W_"+number].ge(.5)
            piece = metadata.loc[keep].copy(); piece["path"] = "Path" + number; parts.append(piece)
        composition(pd.concat(parts, ignore_index=True), "path", f"Monocle2_{orientation}_paths",
                    f"CNV composition by {orientation} root-to-leaf path")

    fig, axes = plt.subplots(1, 3, figsize=(18, 4.8), layout="constrained", gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
    for ax, data, orientation in zip(axes[:2], [main_cells, alt_cells], ["primary", "alternative"]):
        numeric(fig, ax, data.PT_1, f"{orientation}: requested {roots[orientation]['requested_root_group']}", "Pseudotime (a.u.)", domain=pt_domain)
    axes[2].scatter(main_cells.PT_1, alt_cells.PT_1, s=.6, marker="o", linewidths=0,
                    c=metadata.CNV_relative_band.map(COLORS), rasterized=True)
    same = roots["primary"]["actual_root_cell"] == roots["alternative"]["actual_root_cell"]
    rho = float(spearmanr(main_cells.PT_1, alt_cells.PT_1).statistic)
    axes[2].set(xlabel="Primary pseudotime (a.u.)", ylabel="Alternative pseudotime (a.u.)",
                title=f"Same cells | Spearman rho={rho:.3f}" + ("\nSame actual root" if same else ""))
    axes[2].set_box_aspect(span[1]/span[0])
    pd.DataFrame([dict(primary_actual_root=roots["primary"]["actual_root_cell"],
        alternative_actual_root=roots["alternative"]["actual_root_cell"], same_actual_root=same,
        cell_level_descriptive_spearman=rho)]).to_csv(TAB / "root_orientation_comparison.csv", index=False)
    save(fig, "10_root_orientation", "改变计算参考根后，排序是否改变？",
         f"主/替代根是否相同：{same}；全细胞拟时序的描述性Spearman相关为{rho:.4f}。" +
         ("实际根相同，没有形成有效的起点敏感性对照。" if same else "两种根共享DDRTree坐标和树；方向比较不证明真实时间。"),
         ["02_表格/CRC/root_orientation_comparison.csv", "02_表格/CRC/Monocle2_roots.json"])

    cohorts = sorted(metadata.dataset.unique())
    fig, axes = plt.subplots(int(np.ceil(len(cohorts)/4)), 4, figsize=(18, 7.5), layout="constrained")
    for ax, cohort in zip(axes.flat, cohorts):
        keep = metadata.dataset.eq(cohort)
        bands(ax, keep, f"{cohort} | n={keep.sum():,}", legend=False)
    for ax in list(axes.flat)[len(cohorts):]: ax.axis("off")
    fig.legend(handles=[Line2D([], [], marker="o", linestyle="", color=COLORS[b],
        label=f"{b} CNV", markersize=5) for b in BANDS], loc="outside lower center", ncol=3, frameon=False)
    save(fig, "11_cohort_CNV_facets", "共同轨迹中的分支是否由特定队列占据？",
         f"全部{len(cohorts)}队列按同一坐标与树显示，非当前队列为灰色背景；队列分面未重新拟合。", [source_cells])
    metadata.groupby(["dataset", "specimen_id", "sample_id", "CNV_relative_band"], observed=True).size().rename("cells").to_csv(TAB / "frozen_scope_by_unit_CNV.csv")
    pd.read_csv(ROOT / "work/CRC/genes.csv").to_csv(TAB / "ordering_gene_disposition.csv", index=False)

    if args.slingshot:
        u = metadata[["state_UMAP_1", "state_UMAP_2"]].to_numpy()
        for orientation in ["primary", "alternative"]:
            data = aligned(f"Slingshot_{orientation}_cells.csv.gz", metadata)
            curves = pd.read_csv(TAB / f"Slingshot_{orientation}_curves.csv")
            fig, ax = plt.subplots(figsize=(10, 8), layout="constrained")
            ax.scatter(*u[permutation].T, c=metadata.CNV_relative_band.map(COLORS).to_numpy()[permutation],
                       s=.65, marker="o", linewidths=0, rasterized=True)
            for number, curve in curves.groupby("path"):
                ax.plot(curve.UMAP_1, curve.UMAP_2, lw=1, label=f"Path {number}")
            ax.set(title=f"Slingshot {orientation} | fitted in 30D Harmony, projected to UMAP",
                   xlabel="Inherited UMAP 1", ylabel="Inherited UMAP 2")
            ax.set_aspect("equal", adjustable="box")
            ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
            finite = data.filter(regex="^PT_").notna().any(axis=1)
            save(fig, f"12_Slingshot_{orientation}_UMAP", "同细胞Slingshot对照的曲线和CNV分布如何？",
                 f"{finite.sum():,}细胞有至少一条有限拟时序，{(~finite).sum():,}细胞未分配；全部细胞保留显示。曲线为embedCurves投影，二维交叉不证明生物学转换。",
                 [f"02_表格/CRC/Slingshot_{orientation}_cells.csv.gz", f"02_表格/CRC/Slingshot_{orientation}_curves.csv"])

    (LOG / "main_figure_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    (LOG / "main_figure_geometry_checks.json").write_text(json.dumps(dict(
        cells=150917, randomly_omitted_cells=0, identical_axes_in_Monocle2_facets=True,
        full_coordinate_limits=limits, main_figures=len(manifest), visual_acceptance="pending"), indent=2)+"\n")
    print(f"Saved {len(manifest)} PNG/PDF figure pairs; visual inspection still required", flush=True)


if __name__ == "__main__":
    main()
