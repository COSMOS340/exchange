# Reference fragments only. Adapt globals, paths, scope and sample support before use.
# Extracted from the previous Step14 pipeline, not executed by this handoff.


SETS = [
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_MYC_TARGETS_V2",
]

def count_chunks(folder, chunk=2000):
    with h5py.File(folder / "counts.h5") as h:
        shape = h["shape"][:]
        pointers = h["indptr"][:]
        for lo in range(0, int(shape[1]), chunk):
            hi = min(lo+chunk, int(shape[1]))
            a, b = int(pointers[lo]), int(pointers[hi])
            x = sp.csr_matrix((h["data"][a:b].astype(float), h["indices"][a:b],
                              pointers[lo:hi+1]-a), shape=(hi-lo, int(shape[0])))
            yield lo, hi, x

def score_functions():
    members = pd.read_csv(S13 / "输入/Hallmark_members.csv").drop_duplicates()
    assert set(SETS).issubset(set(members.gs_name))
    cfg = json.loads((IN / "analysis_request.json").read_text())
    scope = pd.read_csv(IN / "all_source_cells.csv.gz", low_memory=False)
    eligible = scope[scope.trajectory_eligible].copy()
    eligible["CNV_cohort_percentile"] = eligible.groupby("dataset").cnv_abs_log2.rank(pct=True)
    rank_lookup = eligible.set_index("cell_id").CNV_cohort_percentile
    coverage, validations = [], []
    for spec in cfg["blocks"]:
        key = spec["block"]
        folder = IN / key
        m = pd.read_csv(folder / "metadata.csv.gz", low_memory=False)
        genes = pd.Index(pd.read_csv(folder / "genes.csv").gene)
        panel = sorted(set(members[members.gs_name.isin(SETS)].gene_symbol) & set(genes))
        indices = genes.get_indexer(panel)
        sums, squares = np.zeros(len(panel)), np.zeros(len(panel))
        for lo, hi, x in count_chunks(folder):
            assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(), m.common_gene_UMI.iloc[lo:hi])
            y = x[:, indices].multiply((1e4/m.nCount_input.iloc[lo:hi].to_numpy())[:, None]).tocsr()
            y.data = np.log1p(y.data)
            sums += np.asarray(y.sum(axis=0)).ravel()
            squares += np.asarray(y.power(2).sum(axis=0)).ravel()
        mu = sums/len(m)
        sd = np.sqrt(np.maximum(0, squares/len(m)-mu**2))
        varying = sd > 1e-8
        weights = np.zeros((len(panel), len(SETS)))
        for j, name in enumerate(SETS):
            all_members = set(members.loc[members.gs_name.eq(name), "gene_symbol"])
            hit = np.array([g in all_members for g in panel]) & varying
            assert hit.sum() >= 10 and hit.sum()/len(all_members) >= .5
            weights[hit, j] = 1 / (hit.sum()*sd[hit])
            coverage.append(dict(block=key, pathway=name, resource_genes=len(all_members),
                scored_genes=int(hit.sum()), coverage=float(hit.sum()/len(all_members))))
        shifts = mu @ weights
        scores = np.empty((len(m), len(SETS)), dtype=np.float64)
        check_error = 0.0
        for lo, hi, x in count_chunks(folder):
            y = x[:, indices].multiply((1e4/m.nCount_input.iloc[lo:hi].to_numpy())[:, None]).tocsr()
            y.data = np.log1p(y.data)
            scores[lo:hi] = np.asarray(y @ weights) - shifts
            if lo == 0:
                # A separate dense reconstruction checks scoring on fixed rows.
                dense = x[:30, indices].toarray()
                dense = np.log1p(dense * (1e4/m.nCount_input.iloc[:len(dense)].to_numpy())[:, None])
                for j, name in enumerate(SETS):
                    hit = weights[:, j] > 0
                    expected = ((dense[:, hit]-mu[hit])/sd[hit]).mean(axis=1)
                    check_error = max(check_error, float(np.max(np.abs(expected-scores[:len(dense), j]))))
        assert check_error < 1e-8 and np.isfinite(scores).all()
        d = pd.DataFrame(scores, columns=SETS)
        d.insert(0, "cell_id", m.cell_id.to_numpy())
        d["cnv_abs_log2"] = m.cnv_abs_log2
        d["CNV_cohort_percentile"] = m.cell_id.map(rank_lookup)
        for band in ("Low", "Intermediate", "High"):
            d["CNV_"+band+"_fraction"] = m.CNV_relative_band.eq(band).astype(float)
        assert m.Phase.isin(["G1", "S", "G2M"]).all()
        d["cycling_fraction"] = m.Phase.ne("G1").astype(float)
        d["log1p_UMI"] = np.log1p(m.nCount_input)
        d.to_csv(TAB / key / "function_and_CNV_values.csv.gz", index=False)
        pd.DataFrame(dict(gene=panel, mean=mu, sd=sd, varying=varying)).to_csv(
            TAB / key / "function_gene_standardization.csv", index=False)
        validations.append(dict(block=key, cells=len(m), fixed_row_independent_error=check_error))
        status("function_scores_complete", block=key, cells=len(m))
    pd.DataFrame(coverage).to_csv(TAB / "04_function_gene_coverage.csv", index=False)
    jsave(LOG / "function_score_validation.json", dict(passed=True, checks=validations,
        all_cells_scored=True, denominator="native all-measured-gene UMI",
        interpretation="Mean gene-Z log1p(CP10K), standardized across the joint epithelial pool; not GSEA or flux"))

R_TOPOLOGY = r'''
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(slingshot);library(igraph)})
args=commandArgs(TRUE);input=args[1];out=args[2];key=basename(input)
m=as.data.frame(fread(file.path(input,"metadata.csv.gz")))
h=as.matrix(fread(file.path(input,"Harmony.csv")));p=as.matrix(fread(file.path(input,"PCA.csv")))
rownames(h)=rownames(p)=m$cell_id
edges=function(x,labels){
 s=getLineages(x,clusterLabels=labels)
 e=igraph::as_data_frame(slingMST(s))
 if(nrow(e)==0)return(character())
 sort(apply(e[,c("from","to"),drop=FALSE],1,function(z)paste(sort(z),collapse="|")))
}
main=edges(h,m$coarse_state)
records=list()
compare=function(label,e,detail){
 union_n=length(union(main,e))
 data.frame(block=key,comparison=label,main_edges=length(main),comparison_edges=length(e),
  shared_edges=length(intersect(main,e)),edge_Jaccard=if(union_n)length(intersect(main,e))/union_n else NA_real_,
  detail=detail)
}
records[[1]]=compare("PCA_vs_Harmony",edges(p,m$coarse_state),"Existing full-cell representations; MST only; curves not refitted")
for(ds in sort(unique(m$dataset))){
 keep=m$dataset!=ds
 if(sum(keep)<10)next
 e=tryCatch(edges(h[keep,,drop=FALSE],m$coarse_state[keep]),error=function(err)err)
 if(inherits(e,"error")){
  records[[length(records)+1]]=data.frame(block=key,comparison=paste0("leave_out_",ds),
    main_edges=length(main),comparison_edges=NA,shared_edges=NA,edge_Jaccard=NA,detail=conditionMessage(e))
 }else records[[length(records)+1]]=compare(paste0("leave_out_",ds),e,"MST-only sensitivity; no Harmony or curve refit")
}
fwrite(rbindlist(records),file.path(out,"topology_sensitivity.csv"))
'''

def path_table(key, method, root):
    suffix = "_path_cells" if method == "Monocle2" else "_cells"
    return pd.read_csv(TAB / key / f"{method}_{root}{suffix}.csv.gz")

def dynamics():
    cfg = json.loads((IN / "analysis_request.json").read_text())
    features = ["cnv_abs_log2", "CNV_cohort_percentile",
                "CNV_Low_fraction", "CNV_Intermediate_fraction", "CNV_High_fraction",
                "cycling_fraction", "log1p_UMI", *SETS]
    manifest, endpoints, correlations = [], [], []
    for spec in cfg["blocks"]:
        key = spec["block"]
        m = pd.read_csv(IN / key / "metadata.csv.gz", low_memory=False)
        values = pd.read_csv(TAB / key / "function_and_CNV_values.csv.gz")
        assert np.array_equal(m.cell_id, values.cell_id)
        base = m[["cell_id", "sample_id", "specimen_id", "subject_key", "dataset", "tissue_group"]].copy()
        for feature in features:
            base[feature] = values[feature].to_numpy()
        all_subjects, all_summaries, path_stats = [], [], []
        path_cache = {}
        for method in ("Monocle2", "Slingshot"):
            for root in ("primary", "alternative"):
                p = path_table(key, method, root)
                assert np.array_equal(p.cell_id, m.cell_id)
                path_cache[method, root] = p
                for column in [c for c in p if c.startswith("PT_")]:
                    number = column.split("_")[1]
                    pt = p[column].to_numpy()
                    weight = p["W_"+number].to_numpy()
                    valid = np.isfinite(pt) & (weight >= .5)
                    if not valid.any() or np.ptp(pt[valid]) <= 0:
                        path_stats.append(dict(block=key, method=method, root=root, path=number,
                            cells=int(valid.sum()), status="No varying finite supported pseudotime"))
                        continue
                    scaled = (pt[valid]-pt[valid].min())/np.ptp(pt[valid])
                    d = base.loc[valid].copy()
                    d["bin"] = np.minimum((scaled*12).astype(int), 11)
                    group = ["dataset", "tissue_group", "subject_key", "specimen_id", "sample_id", "bin"]
                    lib = d.groupby(group)[features].mean()
                    lib["cells"] = d.groupby(group).size()
                    lib = lib.reset_index()
                    lib = lib[lib.cells.ge(10)]
                    specimen = lib.groupby(group[:-2]+["bin"])[features].mean().reset_index()
                    subject = specimen.groupby(["dataset", "tissue_group", "subject_key", "bin"])[features].mean().reset_index()
                    long = subject.melt(id_vars=["dataset", "tissue_group", "subject_key", "bin"],
                                        var_name="feature", value_name="value")
                    for name, val in (("block", key), ("method", method), ("root", root), ("path", number)):
                        long[name] = val
                    all_subjects.append(long)
                    summary = long.groupby(["tissue_group", "bin", "feature"]).value.agg(
                        mean="mean", median="median", q25=lambda x:x.quantile(.25),
                        q75=lambda x:x.quantile(.75), subjects="size").reset_index()
                    summary["supported"] = summary.subjects.ge(3)
                    for name, val in (("block", key), ("method", method), ("root", root), ("path", number)):
                        summary[name] = val
                    all_summaries.append(summary)
                    path_stats.append(dict(block=key, method=method, root=root, path=number,
                        cells=int(valid.sum()), subjects=d.subject_key.nunique(),
                        minimum_pt=float(pt[valid].min()), maximum_pt=float(pt[valid].max()),
                        supported_tissue_bins=int(summary[summary.feature.eq("cnv_abs_log2")].supported.sum()),
                        status="Descriptive; shared trunks can occur in multiple paths"))
                    for (tissue, feature), q in summary[summary.supported].groupby(["tissue_group", "feature"]):
                        q = q.sort_values("bin")
                        if len(q) < 2:
                            continue
                        endpoints.append(dict(block=key, method=method, root=root, path=number,
                            tissue=tissue, feature=feature, first_bin=int(q.bin.iloc[0]),
                            last_bin=int(q.bin.iloc[-1]), difference=float(q["mean"].iloc[-1]-q["mean"].iloc[0]),
                            supported_bins=len(q), first_subjects=int(q.subjects.iloc[0]), last_subjects=int(q.subjects.iloc[-1])))
                    for feature in ["cnv_abs_log2", "log1p_UMI", "cycling_fraction"]:
                        y = values.loc[valid, feature].to_numpy()
                        rho = spearmanr(pt[valid], y).statistic if np.ptp(y) > 0 else np.nan
                        correlations.append(dict(block=key, method=method, root=root, path=number,
                            comparator=feature, cells=int(valid.sum()), rho=rho,
                            scope="cell-level descriptive correlation; no P value"))
        for root in ("primary", "alternative"):
            a, b = path_cache["Monocle2", root], path_cache["Slingshot", root]
            for ac in [c for c in a if c.startswith("PT_")]:
                for bc in [c for c in b if c.startswith("PT_")]:
                    valid = (a[ac].notna() & b[bc].notna() &
                             a["W_"+ac.split("_")[1]].ge(.5) & b["W_"+bc.split("_")[1]].ge(.5))
                    if valid.sum() < 30:
                        continue
                    correlations.append(dict(block=key, method="Monocle2_vs_Slingshot", root=root,
                        path=ac+"__"+bc, comparator="pseudotime", cells=int(valid.sum()),
                        rho=spearmanr(a.loc[valid, ac], b.loc[valid, bc]).statistic,
                        scope="same cells and stated directions; not independent validation"))
        pd.concat(all_subjects, ignore_index=True).to_csv(TAB / key / "subject_bin_values.csv.gz", index=False)
        pd.concat(all_summaries, ignore_index=True).to_csv(TAB / key / "path_dynamic_summary.csv", index=False)
        pd.DataFrame(path_stats).to_csv(TAB / key / "path_coverage.csv", index=False)
        manifest.extend(path_stats)
        status("path_dynamics_complete", block=key, paths=len(path_stats))
    pd.DataFrame(manifest).to_csv(TAB / "05_all_path_coverage.csv", index=False)
    pd.DataFrame(endpoints).to_csv(TAB / "06_dynamic_endpoints.csv", index=False)
    pd.DataFrame(correlations).to_csv(TAB / "07_descriptive_correlations.csv", index=False)
    status("all_dynamic_summaries_complete", path_records=len(manifest))

def within_subject_dynamics():
    rows = []
    cfg = json.loads((IN / "analysis_request.json").read_text())
    for spec in cfg["blocks"]:
        key = spec["block"]
        d = pd.read_csv(TAB / key / "subject_bin_values.csv.gz")
        index = ["dataset", "tissue_group", "subject_key", "method", "root", "path", "bin"]
        wide = d.pivot(index=index, columns="feature", values="value").reset_index()
        for group, q in wide.groupby(index[:-1]):
            if len(q) < 4 or q.cnv_abs_log2.nunique() < 3:
                continue
            for feature in SETS:
                if q[feature].nunique() < 3:
                    continue
                rows.append(dict(zip(index[:-1], group), block=key, feature=feature, bins=len(q),
                    rho=float(spearmanr(q.cnv_abs_log2, q[feature]).statistic)))
    result = pd.DataFrame(rows)
    result.to_csv(TAB / "10_within_subject_CNV_function.csv", index=False)
    keys = ["block", "method", "root", "path", "tissue_group", "feature"]
    result.groupby(keys).rho.agg(median="median", q25=lambda x: x.quantile(.25),
        q75=lambda x: x.quantile(.75), subjects="size").reset_index().to_csv(
            TAB / "11_within_subject_CNV_function_summary.csv", index=False)
    status("within_subject_descriptive_cochanges_complete", records=len(result))

def numeric_scatter(ax, xy, values, title, limits=None):
    import matplotlib.pyplot as plt
    v = np.asarray(values, dtype=float)
    finite = np.isfinite(v)
    ax.scatter(xy[:, 0], xy[:, 1], s=.6, color="#dddddd", linewidths=0, rasterized=True)
    if not finite.any():
        clean_axis(ax, title + " | no supported values")
        return
    order = np.flatnonzero(finite)[np.argsort(v[finite], kind="stable")]
    lo, hi = (float(np.min(v[finite])), float(np.max(v[finite]))) if limits is None else limits
    if lo == hi:
        hi = lo + 1
    sc = ax.scatter(xy[order, 0], xy[order, 1], c=v[order], s=1.1, cmap="viridis",
                    vmin=lo, vmax=hi, linewidths=0, rasterized=True)
    clean_axis(ax, title)
    plt.colorbar(sc, ax=ax, shrink=.55, pad=.015)

def category_scatter(ax, xy, values, title, colors=None):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    values = np.asarray(values).astype(str)
    levels = sorted(set(values))
    palette = colors or {level: plt.get_cmap("tab10" if len(levels) <= 10 else "tab20")(i)
                         for i, level in enumerate(levels)}
    for level in levels:
        keep = values == level
        ax.scatter(xy[keep, 0], xy[keep, 1], s=.9, color=palette[level],
                   linewidths=0, rasterized=True)
    clean_axis(ax, title)
    labels = {**TISSUE_NAMES,
              "Prior_mixed_lineage_review": "Prior mixed-lineage review",
              "Known_CD45_enriched_source": "Known CD45 enrichment",
              "All_eligible_cells_fitted_jointly": "All eligible cells fitted jointly"}
    notes = pd.read_csv(S12 / "输入/manual_cluster_identity.csv")
    labels.update({f"C{r.cluster}": f"C{r.cluster} | "+textwrap.fill(r.expression_note, 24)
                   for r in notes.itertuples()})
    handles = [Line2D([], [], marker="o", linestyle="", color=palette[l], markersize=4,
                      label=labels.get(l, l)) for l in levels]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=False, fontsize=7, borderaxespad=0)

def model_geometry(key, method, root="primary"):
    m = pd.read_csv(IN / key / "metadata.csv.gz", low_memory=False)
    if method == "Monocle2":
        d = pd.read_csv(TAB / key / f"Monocle2_{root}_cells.csv.gz")
        assert np.array_equal(m.cell_id, d.cell_id)
        return m, d[["DDRTree_1", "DDRTree_2"]].to_numpy()
    return m, m[["state_UMAP_1", "state_UMAP_2"]].to_numpy()

def add_geometry_lines(ax, key, method, root="primary"):
    if method == "Monocle2":
        centers = pd.read_csv(TAB / key / "Monocle2_centers.csv").set_index("center")
        edges = pd.read_csv(TAB / key / "Monocle2_tree_edges.csv")
        for row in edges.itertuples(index=False, name=None):
            points = centers.loc[[row[0], row[1]], ["DDRTree_1", "DDRTree_2"]].to_numpy()
            ax.plot(points[:, 0], points[:, 1], color="#303030", lw=.6, alpha=.6, zorder=2)
    else:
        curves = pd.read_csv(TAB / key / f"Slingshot_{root}_curves.csv")
        for _, curve in curves.groupby("path"):
            ax.plot(curve.UMAP_1, curve.UMAP_2, color="#303030", lw=.8, alpha=.7, zorder=2)

def overall_main_figures(manifest, methods=("Monocle2", "Slingshot")):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    key = "Overall"
    bands = ["Low", "Intermediate", "High"]
    colors = {"Low": "#3182BD", "Intermediate": "#BDBDBD", "High": "#D94B4B"}
    m = pd.read_csv(IN / key / "metadata.csv.gz", low_memory=False)
    assert len(m) == 223736 and set(m.CNV_relative_band) == set(bands)
    m.groupby(["CNV_relative_band", "tissue_group", "broad_identity"], observed=True).agg(
        cells=("cell_id", "size"), specimens=("specimen_id", "nunique"),
        reported_labels=("subject_key", "nunique")).reset_index().to_csv(
            TAB / "00_overall_CNV_group_composition.csv", index=False)
    geometry_checks = []
    for method in methods:
        _, xy = model_geometry(key, method)
        xspan, yspan = np.ptp(xy, axis=0)
        assert xspan > 0 and yspan > 0
        limits = [(xy[:, j].min()-.04*v, xy[:, j].max()+.04*v)
                  for j, v in enumerate((xspan, yspan))]
        height = max(2.5, min(7.0, 4.8*yspan/xspan))
        handles = [Line2D([], [], color=colors[b], marker="o", linestyle="",
                          markersize=5, label=f"{b} CNV") for b in bands]
        permutation = np.random.default_rng(SEED).permutation(len(m))
        assert np.array_equal(np.sort(permutation), np.arange(len(m)))

        def panel(ax, selected, title, full_color=False):
            if not full_color:
                ax.scatter(xy[:, 0], xy[:, 1], s=.35, color="#ECECEC",
                           linewidths=0, rasterized=True)
            order = permutation[np.asarray(selected)[permutation]]
            point_colors = m.CNV_relative_band.map(colors).to_numpy()[order]
            ax.scatter(xy[order, 0], xy[order, 1], s=.65, c=point_colors,
                       linewidths=0, rasterized=True)
            add_geometry_lines(ax, key, method)
            clean_axis(ax, title)
            ax.set_xlim(limits[0])
            ax.set_ylim(limits[1])
            ax.set_xlabel("DDRTree 1" if method == "Monocle2" else "UMAP 1", fontsize=8)
            ax.set_ylabel("DDRTree 2" if method == "Monocle2" else "UMAP 2", fontsize=8)

        if method == "Monocle2":
            d = pd.read_csv(TAB / key / "Monocle2_primary_cells.csv.gz")
            assert np.array_equal(m.cell_id, d.cell_id) and np.isfinite(d.PT_1).all()
            fig, axes = plt.subplots(1, 2, figsize=(11.5, height+1.4),
                                     gridspec_kw={"wspace": .22})
            panel(axes[0], np.ones(len(m), bool), "A  Relative RNA CNV groups", full_color=True)
            numeric_scatter(axes[1], xy, d.PT_1, "B  Pseudotime | computational reference")
            add_geometry_lines(axes[1], key, method)
            axes[1].set_xlim(limits[0])
            axes[1].set_ylim(limits[1])
            fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
                       bbox_to_anchor=(.5, .01))
            fig.subplots_adjust(bottom=.16, top=.85)
            fig.suptitle(f"Joint epithelial trajectory | all {len(m):,} eligible cells", fontsize=13)
            q = m[["cell_id", "CNV_relative_band", "tissue_group", "subject_key"]].copy()
            q["pseudotime"] = d.PT_1
            stats = q.groupby("CNV_relative_band").pseudotime.agg(
                cells="size", median="median", q25=lambda s: s.quantile(.25),
                q75=lambda s: s.quantile(.75)).reindex(bands)
            stats.to_csv(TAB / "00_overall_CNV_pseudotime_description.csv")
            finding = "；".join(f"{b}组{int(stats.loc[b,'cells']):,}细胞，拟时序中位数"
                                f"{stats.loc[b,'median']:.2f}" for b in bands)
            save_figure(fig, "main/01_joint_Monocle2_CNV_pseudotime", "整体主图：CNV三组与共同拟时序",
                finding+"。两面板来自相同整体模型和坐标；三组未分别拟合。中位数是细胞级描述，"
                "不等于患者级显著差异，也不证明Low向High转换。",
                "02_表格/00_overall_CNV_pseudotime_description.csv", manifest)

        fig, axes = plt.subplots(1, 3, figsize=(15, height+1.2),
                                 gridspec_kw={"wspace": .10})
        counts = []
        for ax, band in zip(axes, bands):
            selected = m.CNV_relative_band.eq(band).to_numpy()
            counts.append(int(selected.sum()))
            panel(ax, selected, f"{band} CNV | {selected.sum():,} cells")
        assert sum(counts) == len(m)
        fig.suptitle(f"{method} | CNV groups on the same joint trajectory", fontsize=13)
        fig.subplots_adjust(top=.84, bottom=.10)
        save_figure(fig, f"main/02_{method}_CNV_facets", f"整体{method}：低、中、高CNV同坐标分面",
            "Low、Intermediate、High分别有"+"/".join(f"{n:,}" for n in counts)+"细胞。"
            "三面板使用相同坐标、范围和完整树线；其他组为浅灰背景。全部细胞参与拟合和显示，没有对子集重拟合。"
            + ("Slingshot线为30维拟合后的UMAP投影，二维交叉不代表生物学转换。" if method == "Slingshot" else ""),
            "02_表格/00_overall_CNV_group_composition.csv", manifest)
        fig, axes = plt.subplots(1, 2, figsize=(11, height+1.4),
                                 gridspec_kw={"wspace": .15})
        for ax, tissue in zip(axes, TISSUES):
            selected = m.tissue_group.eq(tissue).to_numpy()
            panel(ax, selected, f"{TISSUE_NAMES[tissue]} | {selected.sum():,} cells")
        fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
                   bbox_to_anchor=(.5, .01))
        fig.subplots_adjust(top=.85, bottom=.16)
        fig.suptitle(f"{method} | tissue origin on the same joint trajectory", fontsize=13)
        save_figure(fig, f"main/03_{method}_tissue_facets", f"整体{method}：正常/癌旁与CRC分面",
            "正常/癌旁72,819、CRC150,917细胞均来自同一模型；颜色仍为既定CNV三组。"
            "观察来源在分支上的分布，不将组织来源直接当作恶性身份或真实时间。",
            "02_表格/00_overall_CNV_group_composition.csv", manifest)
        geometry_checks.append(dict(method=method, cells=len(m), band_cells=counts,
            equal_axes=True, same_model_all_facets=True, display_subsampled_cells=0,
            x_limits=limits[0], y_limits=limits[1]))
    jsave(LOG / "main_figure_geometry_checks.json", dict(passed=True, checks=geometry_checks))

def validate_dynamics():
    verify_sources()
    cfg = json.loads((IN / "analysis_request.json").read_text())
    checks, outcomes = [], []
    for spec in cfg["blocks"]:
        key = spec["block"]
        assert json.loads((LOG / f"{key}_independent_validation.json").read_text())["passed"]
        m = pd.read_csv(IN / key / "metadata.csv.gz", low_memory=False)
        values = pd.read_csv(TAB / key / "function_and_CNV_values.csv.gz")
        subject = pd.read_csv(TAB / key / "subject_bin_values.csv.gz")
        summary = pd.read_csv(TAB / key / "path_dynamic_summary.csv")
        assert np.array_equal(m.cell_id, values.cell_id)
        features = list(values.columns[1:])
        assert np.isfinite(values[features].to_numpy()).all()
        assert np.max(abs(values[["CNV_Low_fraction", "CNV_Intermediate_fraction",
                                   "CNV_High_fraction"]].sum(axis=1)-1)) == 0
        max_error, tested = 0.0, 0
        for method in ("Monocle2", "Slingshot"):
            complete = json.loads((LOG / f"{key}_{method}_completed.json").read_text())
            assert complete["all_cells_fitted"] and complete["sampled_cells"] == 0 and complete["cells"] == len(m)
            for root in ("primary", "alternative"):
                p = path_table(key, method, root)
                assert np.array_equal(m.cell_id, p.cell_id)
                for col in [c for c in p if c.startswith("PT_")]:
                    number = int(col.split("_")[1])
                    valid = p[col].notna() & p["W_"+str(number)].ge(.5)
                    pt = p.loc[valid, col].to_numpy()
                    if not len(pt) or np.ptp(pt) == 0:
                        continue
                    d = m.loc[valid, ["dataset", "tissue_group", "subject_key", "specimen_id", "sample_id"]].copy()
                    d["bin"] = np.minimum(np.floor(12*(pt-pt.min())/np.ptp(pt)).astype(int), 11)
                    libkeys = ["dataset", "tissue_group", "subject_key", "specimen_id", "sample_id", "bin"]
                    speckeys = ["dataset", "tissue_group", "subject_key", "specimen_id", "bin"]
                    subkeys = ["dataset", "tissue_group", "subject_key", "bin"]
                    cell_n = d.groupby(libkeys).sample_id.transform("size")
                    d = d[cell_n.ge(10)].copy()
                    if d.empty:
                        continue
                    d["inv_cells"] = 1/cell_n.loc[d.index]
                    nlib = d[libkeys].drop_duplicates().groupby(speckeys).size().rename("nlib").reset_index()
                    nspec = d[speckeys].drop_duplicates().groupby(subkeys).size().rename("nspec").reset_index()
                    d["row_index"] = d.index
                    d = d.merge(nlib, on=speckeys, validate="many_to_one").merge(nspec, on=subkeys, validate="many_to_one")
                    weights = d.inv_cells.to_numpy()/d.nlib.to_numpy()/d.nspec.to_numpy()
                    codes, levels = pd.factorize(pd.MultiIndex.from_frame(d[subkeys]), sort=True)
                    aggregate = np.zeros((len(levels), len(features)))
                    np.add.at(aggregate, codes, values.loc[d.row_index, features].to_numpy()*weights[:, None])
                    expected = pd.DataFrame(aggregate, index=levels, columns=features).sort_index()
                    got = subject[subject.method.eq(method) & subject.root.eq(root) & subject.path.eq(number)]
                    actual = got.pivot(index=subkeys, columns="feature", values="value").reindex(
                        index=expected.index, columns=features)
                    assert len(actual) == len(expected) and actual.notna().all().all()
                    error = float(np.max(np.abs(expected.to_numpy()-actual.to_numpy())))
                    assert error < 1e-8
                    max_error = max(max_error, error)
                    tested += expected.size
        regroup = subject.groupby(["method", "root", "path", "tissue_group", "bin", "feature"]).value
        expected = regroup.agg(mean="mean", median="median", subjects="size",
                               q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75))
        actual = summary.set_index(expected.index.names).reindex(expected.index)
        assert np.max(np.abs(expected.to_numpy()-actual[expected.columns].to_numpy())) < 1e-8
        assert np.array_equal(actual.supported, actual.subjects.ge(3))
        groupkeys = ["method", "root", "path", "tissue_group", "bin"]
        fractions = summary[summary.feature.isin(["CNV_Low_fraction", "CNV_Intermediate_fraction",
                                                 "CNV_High_fraction"])].groupby(groupkeys)["mean"].sum()
        assert np.max(np.abs(fractions-1)) < 1e-8
        endpoint = pd.read_csv(TAB / "06_dynamic_endpoints.csv")
        for r in endpoint[endpoint.block.eq(key)].itertuples():
            q = summary[summary.method.eq(r.method) & summary.root.eq(r.root) & summary.path.eq(r.path) &
                        summary.tissue_group.eq(r.tissue) & summary.feature.eq(r.feature) & summary.supported].sort_values("bin")
            assert len(q) == r.supported_bins and q.bin.iloc[0] == r.first_bin and q.bin.iloc[-1] == r.last_bin
            assert abs(q["mean"].iloc[-1]-q["mean"].iloc[0]-r.difference) < 1e-8
        roots = json.loads((TAB / key / "Monocle2_roots.json").read_text())
        obj = pd.read_csv(TAB / key / "Monocle2_objective.csv")
        sling = json.loads((TAB / key / "Slingshot_roots_and_paths.json").read_text())
        for root in ("primary", "alternative"):
            r = roots[root]
            table = pd.read_csv(TAB / key / f"Monocle2_{root}_cells.csv.gz")
            requested_state = r["requested_State"]
            state_cells = int(table.State.astype(str).eq(str(requested_state)).sum()) if isinstance(requested_state, str) else None
            fallback = any("ROOT_FALLBACK:" in line and r["actual_root_cell"] in line
                           for line in (LOG / f"{key}_Monocle2.log").read_text().splitlines())
            outcomes.append(dict(block=key, method="Monocle2", root=root, cells=len(m),
                requested=r["requested_root_group"], actual=r["actual_root_group"], states=r["states"],
                requested_state=requested_state, requested_state_cells=state_cells,
                nearest_leaf_fallback=fallback,
                root_group_matches=r["requested_root_group"] is None or r["requested_root_group"] == r["actual_root_group"],
                objective_iterations=len(obj), hit_iteration_limit=len(obj) >= cfg["monocle_max_iter"],
                final_objective=float(obj.objective.iloc[-1]),
                last_objective_change=float(obj.objective.iloc[-1]-obj.objective.iloc[-2]) if len(obj)>1 else np.nan))
            p = path_table(key, "Slingshot", root)
            pt = [c for c in p if c.startswith("PT_")]
            outcomes.append(dict(block=key, method="Slingshot", root=root, cells=len(m),
                requested=spec["root"] if root == "primary" else spec["alternative"],
                paths=len(pt), unassigned=int(p[pt].isna().all(axis=1).sum()),
                orientation_only=len(m.coarse_state.unique()) == 1 and root == "alternative",
                detail=str(sling[root])))
        checks.append(dict(block=key, all_cells=len(m), independently_checked_values=tested,
                           max_weighted_summary_error=max_error))
        status("independent_dynamics_validated", block=key, values=tested)
    pd.DataFrame(outcomes).to_csv(TAB / "09_model_reference_audit.csv", index=False)
    jsave(LOG / "dynamic_validation.json", dict(passed=True, checks=checks,
        method="Independent cell-level weights reconstruct library/specimen/subject means on every supported path",
        no_cell_level_inference=True))
