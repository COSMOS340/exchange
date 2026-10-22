"""Path dynamics with library/specimen/reported-label weighting and explicit gaps."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/02_表格/AllEpithelial"
UNITS = TAB / "dynamics"
FIG = ROOT / "results/01_图形/dynamics"
LOG = ROOT / "results/04_日志和检查"
SETS = json.loads((ROOT / "config.json").read_text())["function_sets"]
FEATURES = ["cnv_abs_log2", "CNV_Low_fraction", "CNV_Intermediate_fraction", "CNV_High_fraction",
            *SETS, "cycling_fraction", "log1p_UMI"]


def load_aligned(path, metadata):
    values = pd.read_csv(path)
    assert values.cell_id.is_unique and set(values.cell_id) == set(metadata.cell_id)
    return values.set_index("cell_id").loc[metadata.cell_id].reset_index()


def plot_path(summary, unsupported, key, coverage):
    result = []
    for proxy in [False, True]:
        if proxy and unsupported.empty:
            continue
        fig, axes = plt.subplots(4, 3, figsize=(15, 13), layout="constrained")
        x = (np.arange(12) + .5) / 12
        for ax, feature in zip(axes.flat, FEATURES):
            if proxy:
                for _, group in unsupported.groupby(["dataset", "specimen_id"]):
                    values = group.set_index("bin")[feature].reindex(range(12))
                    ax.plot(x, values, marker="o", markersize=3, linewidth=.8, alpha=.65)
                finite = unsupported[feature].notna().any()
            else:
                group = summary[summary.feature.eq(feature)].set_index("bin").reindex(range(12))
                ax.plot(x, group["mean"], marker="o", markersize=3, linewidth=1.4, color="#2C628A")
                ax.fill_between(x, group.q25, group.q75, color="#2C628A", alpha=.2)
                finite = group["mean"].notna().any()
            title = feature.removeprefix("HALLMARK_").replace("_", " ")
            ax.set(title=title, xlabel="Within-path scaled pseudotime (0-1)", xlim=(0, 1))
            ax.tick_params(labelsize=7)
            ax.spines[["top", "right"]].set_visible(False)
            if "fraction" in feature:
                ax.set(ylim=(-.03, 1.03), ylabel="Proportion (0-1)")
            elif feature in SETS:
                ax.set_ylabel("Mean gene Z of log1p(CP10K)")
            elif feature == "cnv_abs_log2":
                ax.set_ylabel("Mean absolute log2 relative RNA CNV")
            else:
                ax.set_ylabel("log1p(native total UMI)")
            if not finite:
                ax.text(.5, .5, "Insufficient supported coverage", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8)
        label = "Unsupported patient links: individual specimen means" if proxy else "Supported reported labels: equal-label mean and IQR"
        fig.suptitle(f"{key} | {label}\nEach library/bin needs >=10 cells; gaps are retained; descriptive only", fontsize=12)
        suffix = "proxy_specimens" if proxy else "supported_labels"
        name = f"{key}_{suffix}"
        fig.canvas.draw()
        fig.set_layout_engine("none")
        fig.savefig(FIG / (name + ".png"), dpi=160, bbox_inches="tight")
        fig.savefig(FIG / (name + ".pdf"), dpi=160, bbox_inches="tight")
        plt.close(fig)
        cnv = summary[summary.feature.eq("cnv_abs_log2")]
        supported_bins = int(cnv.supported.sum())
        finding = (f"该路径纳入{coverage['cells']:,}个符合路径成员/权重规则的细胞；"
                   f"连续CNV有{supported_bins}/12箱满足至少3个有支持报告标签。"
                   "空箱保留，曲线不跨空箱连接；不同箱的标签组成可能不同，不等于同一患者进展。")
        if proxy:
            finding = (f"单独显示{unsupported[['dataset','specimen_id']].drop_duplicates().shape[0]}个缺少患者对应支持的标本；"
                       "每条线是一个标本，不纳入主体报告标签均值，不声称独立患者。")
        result.append(dict(name=name, png=f"01_图形/dynamics/{name}.png", pdf=f"01_图形/dynamics/{name}.pdf",
            question="路径内CNV与表达活动在有覆盖支持的分箱中如何变化？", finding=finding,
            statistical_unit="individual proxy specimens" if proxy else "equal supported reported labels after equal library/specimen weighting",
            sources=[f"02_表格/AllEpithelial/dynamics/{key}_summary.csv",
                     f"02_表格/AllEpithelial/dynamics/{key}_specimen_bins.csv.gz",
                     f"02_表格/AllEpithelial/dynamics/{key}_reported_label_bins.csv.gz"], visual_review="pending"))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--method', choices=['all', 'Monocle2', 'Slingshot'], default='all')
    selected = parser.parse_args().method
    methods = ['Monocle2', 'Slingshot'] if selected == 'all' else [selected]
    for method in methods:
        if (LOG / f"{method}_dynamics_manifest.json").exists():
            raise RuntimeError(f"{method} dynamics exist; validate or version before recomputing")
    UNITS.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(ROOT / "work/AllEpithelial/metadata.csv.gz", low_memory=False)
    assert len(metadata) == 223736 and metadata.cell_id.is_unique
    assert metadata.reported_patient_link_supported.isin([True, False]).all()
    values = load_aligned(TAB / "function_and_CNV_values.csv.gz", metadata)
    keys = ["cell_id", "dataset", "subject_key", "specimen_id", "capture_id", "reported_patient_link_supported"]
    base = metadata[keys].copy()
    base[FEATURES] = values[FEATURES].to_numpy()
    source_tables = {}
    manifest, coverage_records, all_summaries, changes, comparisons = [], [], [], [], []
    for method in methods:
        for orientation in ["primary", "alternative"]:
            suffix = "path_cells" if method == "Monocle2" else "cells"
            paths = load_aligned(TAB / f"{method}_{orientation}_{suffix}.csv.gz", metadata)
            source_tables[method, orientation] = paths
            for pt_name in [name for name in paths if name.startswith("PT_")]:
                number = pt_name.split("_")[1]
                weight = paths["W_" + number].to_numpy()
                pt = paths[pt_name].to_numpy()
                valid = np.isfinite(pt) & (weight >= .5)
                key = f"{method}_{orientation}_Path{number}"
                coverage = dict(method=method, root=orientation, path=number, cells=int(valid.sum()),
                                finite_time_cells=int(np.isfinite(pt).sum()), status="pending")
                if not valid.any() or np.ptp(pt[valid]) <= 0:
                    coverage.update(status="not_assessable_no_varying_supported_pseudotime")
                    coverage_records.append(coverage); continue
                minimum, maximum = float(pt[valid].min()), float(pt[valid].max())
                scaled = (pt[valid]-minimum)/(maximum-minimum)
                bins = np.minimum(np.floor(12*scaled).astype(int), 11)
                assert bins.min() >= 0 and bins.max() <= 11
                frame = base.loc[valid].copy()
                frame["bin"] = bins
                assignments = frame[keys + ["bin"]].copy()
                assignments["raw_pseudotime"] = pt[valid]
                assignments["scaled_pseudotime"] = scaled
                assignments["weight"] = weight[valid]
                assignments.to_csv(UNITS / f"{key}_cell_bins.csv.gz", index=False)
                pd.DataFrame(dict(bin=np.arange(12), scaled_left=np.arange(12)/12,
                    scaled_right=np.arange(1,13)/12, raw_left=minimum+np.arange(12)/12*(maximum-minimum),
                    raw_right=minimum+np.arange(1,13)/12*(maximum-minimum),
                    right_inclusive=[False]*11+[True])).to_csv(UNITS / f"{key}_boundaries.csv", index=False)
                libkeys = ["dataset", "subject_key", "specimen_id", "capture_id", "reported_patient_link_supported", "bin"]
                grouped = frame.groupby(libkeys, dropna=False, observed=True)
                libraries = grouped[FEATURES].mean().join(grouped.size().rename("cells")).reset_index()
                libraries["eligible"] = libraries.cells.ge(10)
                libraries.to_csv(UNITS / f"{key}_library_bins.csv.gz", index=False)
                eligible = libraries[libraries.eligible]
                speckeys = ["dataset", "subject_key", "specimen_id", "reported_patient_link_supported", "bin"]
                grouped = eligible.groupby(speckeys, dropna=False, observed=True)
                specimens = grouped[FEATURES].mean().join(grouped.cells.sum().rename("cells")).join(
                    grouped.size().rename("libraries")).reset_index()
                specimens.to_csv(UNITS / f"{key}_specimen_bins.csv.gz", index=False)
                supported = specimens[specimens.reported_patient_link_supported]
                unsupported = specimens[~specimens.reported_patient_link_supported]
                unsupported.to_csv(UNITS / f"{key}_unsupported_specimens.csv.gz", index=False)
                labelkeys = ["dataset", "subject_key", "bin"]
                grouped = supported.groupby(labelkeys, dropna=False, observed=True)
                labels = grouped[FEATURES].mean().join(grouped.cells.sum().rename("cells")).join(
                    grouped.libraries.sum().rename("libraries")).join(grouped.size().rename("specimens")).reset_index()
                labels.to_csv(UNITS / f"{key}_reported_label_bins.csv.gz", index=False)
                long = labels.melt(id_vars=labelkeys, value_vars=FEATURES, var_name="feature", value_name="value")
                summary = long.groupby(["bin", "feature"], observed=True).value.agg(
                    raw_mean="mean", raw_q25=lambda x: x.quantile(.25), raw_q75=lambda x: x.quantile(.75), reported_labels="count")
                complete = pd.MultiIndex.from_product([range(12), FEATURES], names=["bin", "feature"])
                summary = summary.reindex(complete).reset_index()
                summary["reported_labels"] = summary.reported_labels.fillna(0).astype(int)
                summary["supported"] = summary.reported_labels.ge(3)
                for name in ["mean", "q25", "q75"]:
                    summary[name] = summary["raw_"+name].where(summary.supported)
                for name, value in [("method", method), ("root", orientation), ("path", number)]:
                    summary[name] = value
                summary.to_csv(UNITS / f"{key}_summary.csv", index=False)
                all_summaries.append(summary)
                # Within-label changes preserve whose two observed bins are compared.
                for (dataset, subject), block in labels.groupby(["dataset", "subject_key"]):
                    block = block.sort_values("bin")
                    for feature in FEATURES:
                        finite = block[block[feature].notna()]
                        if len(finite) < 2: continue
                        changes.append(dict(method=method, root=orientation, path=number, dataset=dataset,
                            subject_key=subject, feature=feature, observed_bins=len(finite),
                            first_bin=int(finite.bin.iloc[0]), last_bin=int(finite.bin.iloc[-1]),
                            first_value=float(finite[feature].iloc[0]), last_value=float(finite[feature].iloc[-1]),
                            difference=float(finite[feature].iloc[-1]-finite[feature].iloc[0])))
                coverage.update(status="descriptive", minimum_pt=minimum, maximum_pt=maximum,
                    eligible_library_bins=len(eligible), supported_reported_labels=labels[["dataset","subject_key"]].drop_duplicates().shape[0],
                    unsupported_specimens=unsupported[["dataset","specimen_id"]].drop_duplicates().shape[0],
                    supported_CNV_bins=int(summary[summary.feature.eq("cnv_abs_log2")].supported.sum()))
                coverage_records.append(coverage)
                manifest.extend(plot_path(summary, unsupported, key, coverage))
                print(key, coverage, flush=True)
    for method in methods:
        pd.concat([x for x in all_summaries if x.method.iloc[0] == method], ignore_index=True).to_csv(
            TAB / f'{method}_path_dynamic_summary.csv', index=False)
        pd.DataFrame([x for x in coverage_records if x['method'] == method]).to_csv(
            TAB / f'{method}_path_coverage.csv', index=False)
        pd.DataFrame([x for x in changes if x['method'] == method]).to_csv(
            TAB / f'{method}_within_reported_label_changes.csv.gz', index=False)
        items = [x for x in manifest if x['name'].startswith(method+'_')]
        (LOG / f'{method}_dynamics_manifest.json').write_text(json.dumps(items, ensure_ascii=False, indent=2)+'\n')
    # Model-independent postprocessing may finish for Monocle2 while Slingshot fits.
    # Assemble the combined outputs only once both full-cell methods are available.
    if not all((LOG / f'{method}_dynamics_manifest.json').exists() for method in ['Monocle2','Slingshot']):
        print(f'Saved {len(manifest)} PNG/PDF pairs for {selected}; other method still pending', flush=True)
        return
    for method in ['Monocle2','Slingshot']:
        for orientation in ['primary','alternative']:
            if (method, orientation) not in source_tables:
                suffix = 'path_cells' if method == 'Monocle2' else 'cells'
                source_tables[method, orientation] = load_aligned(
                    TAB / f'{method}_{orientation}_{suffix}.csv.gz', metadata)
    for orientation in ["primary", "alternative"]:
        a, b = source_tables["Monocle2", orientation], source_tables["Slingshot", orientation]
        for ac in [c for c in a if c.startswith("PT_")]:
            for bc in [c for c in b if c.startswith("PT_")]:
                valid = a[ac].notna() & b[bc].notna() & a["W_"+ac.split("_")[1]].ge(.5) & b["W_"+bc.split("_")[1]].ge(.5)
                rho = np.nan
                if valid.sum() >= 30 and a.loc[valid,ac].nunique()>1 and b.loc[valid,bc].nunique()>1:
                    rho = float(spearmanr(a.loc[valid,ac], b.loc[valid,bc]).statistic)
                comparisons.append(dict(root=orientation, monocle_path=ac, slingshot_path=bc,
                    overlapping_supported_cells=int(valid.sum()), descriptive_spearman=rho,
                    direction="stored computational-reference direction; no pair selection", independent_validation=False))
    for suffix in ['path_dynamic_summary.csv','path_coverage.csv','within_reported_label_changes.csv.gz']:
        pd.concat([pd.read_csv(TAB/f'{method}_{suffix}') for method in ['Monocle2','Slingshot']],
                  ignore_index=True).to_csv(TAB/suffix, index=False)
    pd.DataFrame(comparisons).to_csv(TAB / "all_path_method_comparisons.csv", index=False)
    manifest = sum([json.loads((LOG/f'{method}_dynamics_manifest.json').read_text())
                    for method in ['Monocle2','Slingshot']], [])
    (LOG / "dynamics_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n")
    print(f"Saved {len(manifest)} PNG/PDF pairs; labels are reported labels, not verified independent patients", flush=True)


if __name__ == "__main__":
    main()
