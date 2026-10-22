"""Score six frozen gene sets on every prepared epithelial cell using native total UMI."""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/02_表格/AllEpithelial"
LOG = ROOT / "results/04_日志和检查"
SETS = json.loads((ROOT / "config.json").read_text())["function_sets"]


def count_chunks(chunk=2000):
    # The stored matrix is genes x cells CSC; the same arrays represent its
    # cells x genes CSR transpose. Chunks are scoring I/O, never separate fits.
    with h5py.File(ROOT / "work/AllEpithelial/counts.h5") as h:
        shape = h["shape"][:]
        pointers = h["indptr"][:]
        for lo in range(0, int(shape[1]), chunk):
            hi = min(lo + chunk, int(shape[1]))
            a, b = int(pointers[lo]), int(pointers[hi])
            matrix = sp.csr_matrix((h["data"][a:b].astype(float), h["indices"][a:b],
                                   pointers[lo:hi+1] - a), shape=(hi-lo, int(shape[0])))
            yield lo, hi, matrix


def main():
    output = TAB / "function_and_CNV_values.csv.gz"
    if output.exists():
        raise RuntimeError("Scoring output exists; validate it or use a new output version")
    metadata = pd.read_csv(ROOT / "work/AllEpithelial/metadata.csv.gz", low_memory=False)
    assert len(metadata) == 223736 and metadata.cell_id.is_unique
    assert np.isfinite(metadata.nCount_input).all() and metadata.nCount_input.gt(0).all()
    members = pd.read_csv(ROOT / "data/gene_sets/Hallmark_members.csv").drop_duplicates()
    members = members[members.gs_name.isin(SETS)]
    assert set(members.gs_name) == set(SETS)
    genes = pd.Index(pd.read_csv(ROOT / "work/AllEpithelial/genes.csv").gene)
    panel = sorted(set(members.gene_symbol) & set(genes))
    indices = genes.get_indexer(panel)
    sums, squares = np.zeros(len(panel)), np.zeros(len(panel))
    for lo, hi, x in count_chunks():
        assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(), metadata.common_gene_UMI.iloc[lo:hi])
        y = x[:, indices].multiply((1e4 / metadata.nCount_input.iloc[lo:hi].to_numpy())[:, None]).tocsr()
        y.data = np.log1p(y.data)
        sums += np.asarray(y.sum(axis=0)).ravel()
        squares += np.asarray(y.power(2).sum(axis=0)).ravel()
    mean = sums / len(metadata)
    sd = np.sqrt(np.maximum(0, squares / len(metadata) - mean**2))
    variable = sd > 1e-8
    weights = np.zeros((len(panel), len(SETS)))
    coverage, membership_rows = [], []
    valid_sets = []
    panel_index = {gene: index for index, gene in enumerate(panel)}
    for j, name in enumerate(SETS):
        original = set(members.loc[members.gs_name.eq(name), "gene_symbol"])
        hit = np.array([gene in original for gene in panel]) & variable
        valid = bool(hit.sum() >= 10 and hit.sum() / len(original) >= .5)
        valid_sets.append(valid)
        if valid:
            weights[hit, j] = 1 / (int(hit.sum()) * sd[hit])
        coverage.append(dict(pathway=name, resource_genes=len(original),
            variable_hit_genes=int(hit.sum()), coverage=float(hit.sum()/len(original)),
            assessable=valid, status="assessable" if valid else "insufficient_variable_gene_coverage"))
        for gene in sorted(original):
            present = gene in panel_index
            varies = bool(variable[panel_index[gene]]) if present else False
            membership_rows.append(dict(pathway=name, gene=gene, present=present,
                variable=varies, used_for_score=valid and varies,
                status="used" if valid and varies else "missing" if not present else "nonvariable" if not varies else "set_not_assessable"))
    shift = mean @ weights
    scores = np.full((len(metadata), len(SETS)), np.nan)
    independent_error = 0.
    for lo, hi, x in count_chunks():
        y = x[:, indices].multiply((1e4 / metadata.nCount_input.iloc[lo:hi].to_numpy())[:, None]).tocsr()
        y.data = np.log1p(y.data)
        block = np.asarray(y @ weights) - shift
        block[:, ~np.array(valid_sets)] = np.nan
        scores[lo:hi] = block
        if lo == 0:
            # Fixed-row arithmetic QA only; all 150,917 cells were scored.
            dense = np.log1p(x[:30, indices].toarray() * (1e4 / metadata.nCount_input.iloc[:30].to_numpy())[:, None])
            for j, valid in enumerate(valid_sets):
                if valid:
                    hit = weights[:, j] > 0
                    expected = ((dense[:, hit] - mean[hit]) / sd[hit]).mean(axis=1)
                    independent_error = max(independent_error, float(np.max(np.abs(expected-scores[:30, j]))))
    assert independent_error < 1e-8
    if any(valid_sets):
        assert np.isfinite(scores[:, valid_sets]).all()
        assert np.abs(scores[:, valid_sets].mean(axis=0)).max() < 1e-8
    values = pd.DataFrame(scores, columns=SETS)
    values.insert(0, "cell_id", metadata.cell_id.to_numpy())
    values["cnv_abs_log2"] = pd.read_csv(ROOT / "work/AllEpithelial/metadata.csv.gz", usecols=["cnv_abs_log2"], dtype=str)["cnv_abs_log2"]
    for band in ["Low", "Intermediate", "High"]:
        values["CNV_" + band + "_fraction"] = metadata.CNV_relative_band.eq(band).astype(float)
    assert metadata.Phase.isin(["G1", "S", "G2M"]).all()
    values["cycling_fraction"] = metadata.Phase.ne("G1").astype(float)
    values["log1p_UMI"] = np.log1p(metadata.nCount_input)
    values.to_csv(output, index=False)
    pd.DataFrame(dict(gene=panel, mean_log1p_CP10K=mean, population_sd=sd,
                      variable=variable)).to_csv(TAB / "function_gene_standardization.csv", index=False)
    pd.DataFrame(coverage).to_csv(TAB / "function_gene_coverage.csv", index=False)
    pd.DataFrame(membership_rows).to_csv(TAB / "function_gene_membership.csv", index=False)
    result = dict(passed=True, cells=len(metadata), all_cells_scored=True,
        denominator="nCount_input: original all-measured-gene UMI", pool="all 223736 epithelial cells",
        population_standard_deviation=True, independently_reconstructed_rows=30,
        independent_max_error=independent_error, coverage=coverage,
        interpretation="Mean gene Z of log1p(CP10K); expression activity, not flux or new GSEA")
    (LOG / "function_score_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
