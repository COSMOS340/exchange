"""Portable preparation and fitting for the CRC-only Monocle2 handoff.

No fit runs without the explicit fit command. The scientific deliverable also
requires the diagnostics, figures and human review specified in TASK_SPEC.md.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
EXPECTED_SOURCE = 223736
EXPECTED_CRC = 150917


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024**2), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def manifest(folder):
    return [{"path": p.relative_to(folder).as_posix(), "bytes": p.stat().st_size,
             "sha256": digest(p)} for p in sorted(folder.rglob("*")) if p.is_file()]


def verify(folder, records):
    for item in records:
        p = (folder / item["path"]).resolve()
        if not p.is_relative_to(folder.resolve()):
            raise ValueError("Unsafe manifest path")
        if not p.is_file() or p.stat().st_size != item["bytes"] or digest(p) != item["sha256"]:
            raise ValueError(f"File verification failed: {item['path']}")


def crc_mask(metadata):
    import numpy as np
    required = ["cell_id", "tissue_group", "trajectory_eligible", "cnv_summary_eligible",
                "identity_review", "known_CD45_enriched", "CNV_relative_band"]
    if not set(required).issubset(metadata.columns):
        raise ValueError("Required input fields are missing")
    if not metadata.cell_id.is_unique or metadata.cell_id.isna().any():
        raise ValueError("Invalid cell IDs")
    for name in ["trajectory_eligible", "cnv_summary_eligible", "identity_review", "known_CD45_enriched"]:
        if metadata[name].isna().any() or not metadata[name].isin([True, False]).all():
            raise ValueError(f"Invalid boolean field: {name}")
    keep = metadata.tissue_group.eq("Primary_CRC")
    selected = metadata.loc[keep]
    if not selected.trajectory_eligible.all() or not selected.cnv_summary_eligible.all():
        raise ValueError("CRC selection contains ineligible cells")
    if selected.identity_review.any() or selected.known_CD45_enriched.any():
        raise ValueError("CRC selection violates the frozen source rules")
    if not selected.CNV_relative_band.isin(["Low", "Intermediate", "High"]).all():
        raise ValueError("Unrecognized CNV band; do not guess or relabel")
    return np.asarray(keep)


def choose_roots(metadata, harmony):
    import numpy as np
    if len(metadata) != len(harmony) or not np.isfinite(harmony).all():
        raise ValueError("Invalid Harmony coordinates")
    coverage = metadata.groupby(["coarse_state", "specimen_id"]).size().rename("cells").reset_index()
    support = coverage[coverage.cells.ge(30)].groupby("coarse_state").specimen_id.nunique()
    counts = metadata.coarse_state.value_counts()
    labels = sorted(counts.index, key=lambda s: (-support.get(s, 0), -counts[s], s))
    if len(labels) < 2:
        raise ValueError("At least two inherited clusters are required for this root comparison")
    primary = labels[0]
    centers = {s: harmony[metadata.coarse_state.eq(s).to_numpy()].mean(axis=0) for s in labels}
    pool = [s for s in labels if s != primary and support.get(s, 0) >= 5]
    pool = pool or [s for s in labels if s != primary]
    alternative = max(pool, key=lambda s: float(np.sum((centers[s] - centers[primary])**2)))
    return str(primary), str(alternative), coverage


def build_engine(source):
    start = source.index("patch_monocle2_dplyr <- function()")
    engine = source[start:]
    old = "patch_sparse_size_factors();patch_blocked_dispersion()"
    if engine.count(old) != 1:
        raise ValueError("Unexpected source engine; do not guess a replacement")
    engine = engine.replace(";" + old, "")
    if "patch_blocked_dispersion" in engine or "patch_sparse_size_factors" in engine:
        raise ValueError("Mac blocked preprocessing was not removed")
    return engine


def verify_data():
    records = json.loads((ROOT / "INPUT_MANIFEST.json").read_text(encoding="utf-8"))
    verify(ROOT / "data", records)
    print(f"Verified {len(records)} input files", flush=True)


def prepare():
    import h5py
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    verify_data()
    work = ROOT / "work"
    if work.exists():
        raise RuntimeError("work already exists. Verify it or use a new version; never overwrite it blindly.")
    source = ROOT / "data/Overall"
    m = pd.read_csv(source / "metadata.csv.gz", low_memory=False)
    if len(m) != EXPECTED_SOURCE:
        raise ValueError("Unexpected source cell count")
    keep = crc_mask(m)
    if int(keep.sum()) != EXPECTED_CRC:
        raise ValueError("Unexpected CRC cell count")
    selected = m.loc[keep].copy().reset_index(drop=True)
    genes = pd.read_csv(source / "genes.csv")
    if not genes.gene.is_unique or genes.gene.isna().any():
        raise ValueError("Gene identifiers must be unique and present")
    if genes.ordering.isna().any() or not genes.ordering.isin([True, False]).all():
        raise ValueError("Invalid ordering-gene flags")
    with h5py.File(source / "counts.h5", "r") as h:
        shape = tuple(map(int, h["shape"][:]))
        if shape != (len(genes), len(m)):
            raise ValueError("Matrix dimensions do not match the metadata")
        x = sp.csc_matrix((h["data"][:], h["indices"][:], h["indptr"][:]), shape=shape)
    if not x.has_canonical_format or (x.data < 0).any():
        raise ValueError("Input counts are not canonical nonnegative counts")
    if not np.issubdtype(x.data.dtype, np.integer):
        raise ValueError("Expected integer raw counts")
    x = x[:, keep].tocsc()
    if not np.array_equal(np.asarray(x.sum(axis=0)).ravel(), selected.common_gene_UMI.to_numpy()):
        raise ValueError("Cell order or common-gene UMI totals differ")
    detected = np.asarray((x > 0).sum(axis=1)).ravel()
    original_ordering = genes.ordering.to_numpy(dtype=bool)
    genes["source_ordering"] = original_ordering
    genes["detected_cells"] = detected
    genes["ordering"] = original_ordering & (detected >= 10)
    if genes.ordering.sum() < 2:
        raise ValueError("Insufficient ordering genes")
    harmony = np.load(source / "Harmony.npy", allow_pickle=False)
    pca = np.load(source / "PCA.npy", allow_pickle=False)
    if harmony.shape != (len(m), 30) or pca.shape != (len(m), 30) or not np.isfinite(pca).all():
        raise ValueError("Unexpected embedding dimensions")
    harmony, pca = harmony[keep], pca[keep]
    primary, alternative, coverage = choose_roots(selected, harmony)
    selected["block"] = "CRC"
    selected["trajectory_disposition"] = "All_eligible_CRC_cells_fitted_jointly"
    target = work / "CRC"
    target.mkdir(parents=True)
    selected.to_csv(target / "metadata.csv.gz", index=False)
    genes.to_csv(target / "genes.csv", index=False)
    coverage.to_csv(work / "root_coverage.csv", index=False)
    with h5py.File(target / "counts.h5", "w") as h:
        for name, value in [("data", x.data), ("indices", x.indices), ("indptr", x.indptr)]:
            h.create_dataset(name, data=value, compression="gzip")
        h.create_dataset("shape", data=np.asarray(x.shape, dtype=np.int64))
    pd.DataFrame(harmony).to_csv(target / "Harmony.csv", index=False)
    pd.DataFrame(pca).to_csv(target / "PCA.csv", index=False)
    cfg = json.loads((ROOT / "config.json").read_text(encoding="utf-8"))
    cfg["blocks"] = [{"block": "CRC", "identity": "CRC-origin epithelial cells",
                      "cells": len(selected), "root": primary, "alternative": alternative}]
    save_json(work / "analysis_request.json", cfg)
    source_code = (ROOT / "reference/Overall_Monocle2_generated.R").read_text(encoding="utf-8")
    (work / "fit_generated.R").write_text(build_engine(source_code), encoding="utf-8")
    summary = {"source_cells": len(m), "CRC_cells": len(selected),
               "normal_cells_not_fitted": int((~keep).sum()), "sampled_cells": 0,
               "genes": len(genes), "ordering_genes": int(genes.ordering.sum()),
               "cohorts": selected.dataset.nunique(), "specimens": selected.specimen_id.nunique(),
               "libraries": selected.sample_id.nunique(), "reported_subject_labels": selected.subject_key.nunique(),
               "primary_reference_cluster": primary, "alternative_reference_cluster": alternative,
               "CNV_bands_unchanged": True, "UMI_and_cell_order_verified": True,
               "preparation_only": True}
    save_json(work / "preparation_summary.json", summary)
    save_json(ROOT / "WORK_MANIFEST.json", manifest(work))
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


def fit(rscript, method):
    verify_data()
    records = json.loads((ROOT / "WORK_MANIFEST.json").read_text(encoding="utf-8"))
    verify(ROOT / "work", records)
    executable = shutil.which(rscript) if not Path(rscript).is_file() else str(Path(rscript).resolve())
    if executable is None:
        raise FileNotFoundError("Rscript not found; supply its verified path with --rscript")
    env = os.environ.copy()
    env.pop("R_MAX_VSIZE", None)
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        env[key] = "2"
    preflight = """stopifnot(getRversion()>='4.5.0');
pkgs=c('monocle','DDRTree','igraph','Rcpp','Matrix','data.table','jsonlite','rhdf5','dplyr','rlang','R.utils');
stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)));
stopifnot(as.character(packageVersion('monocle'))=='2.38.0',
as.character(packageVersion('DDRTree'))=='0.1.6',
as.character(packageVersion('igraph'))=='2.2.2');
Rcpp::cppFunction('int compiler_check(){return 1;}');stopifnot(compiler_check()==1);sessionInfo()"""
    if method == "Slingshot":
        preflight += ";stopifnot(requireNamespace('slingshot',quietly=TRUE))"
    result = subprocess.run([executable, "-e", preflight], env=env, text=True, capture_output=True)
    print(result.stdout, end="")
    if result.returncode:
        raise RuntimeError(result.stderr + "\nEnvironment differs; establish an isolated compatible environment before fitting.")
    output = ROOT / "results"
    obj, tab, log = output / "03_对象/CRC", output / "02_表格/CRC", output / "04_日志和检查"
    for p in [obj, tab, log]:
        p.mkdir(parents=True, exist_ok=True)
    binding = {"method": method, "work_manifest_sha256": digest(ROOT / "WORK_MANIFEST.json"),
               "Rscript": executable, "code_sha256": digest(ROOT / "work/fit_generated.R"),
               "preflight_sha256": hashlib.sha256(result.stdout.encode()).hexdigest()}
    binding_file = log / f"{method}_binding.json"
    if binding_file.exists() and json.loads(binding_file.read_text(encoding="utf-8")) != binding:
        raise RuntimeError("Checkpoint binding differs; do not reuse these outputs")
    if not binding_file.exists() and (obj / "Monocle2_DDRTree_checkpoint.rds").exists() and method == "Monocle2":
        raise RuntimeError("Unbound checkpoint found")
    save_json(binding_file, binding)
    (log / f"{method}_preflight.txt").write_text(result.stdout + result.stderr, encoding="utf-8")
    completed = log / f"CRC_{method}_completed.json"
    if completed.exists():
        raise RuntimeError("A completion marker exists; verify it rather than silently rerunning")
    command = [executable, str(ROOT / "work/fit_generated.R"), str(ROOT / "work/CRC"),
               str(obj), str(tab), str(log), str(ROOT / "work/analysis_request.json"), method]
    print("Starting full-cell " + method + "; inspect results/04_日志和检查 for progress.", flush=True)
    with open(log / f"{method}_console.log", "a", encoding="utf-8") as f:
        result = subprocess.run(command, env=env, stdout=f, stderr=subprocess.STDOUT)
    save_json(log / f"{method}_exit.json", {"returncode": result.returncode, "method": method})
    if result.returncode or not completed.exists():
        raise RuntimeError("Fit did not complete. Preserve logs; do not change the cell scope.")
    print("Fit saved. Numeric checks, figures and scientific review are still required.", flush=True)


def pack_return(label):
    if not label or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in label):
        raise ValueError("Use a simple new ASCII version label")
    results = ROOT / "results"
    if not results.exists() or not any(results.rglob("*")):
        raise ValueError("No results to return")
    summary = results / "RUN_SUMMARY.md"
    if not summary.exists():
        raise ValueError("Write results/RUN_SUMMARY.md with completed/failed/pending items before returning")
    target = ROOT / "return" / label
    target.mkdir(parents=True, exist_ok=False)
    save_json(target / "RESULT_MANIFEST.json", manifest(results))
    archive = target / "results.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1, allowZip64=True) as z:
        for p in sorted(results.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(ROOT).as_posix())
        for p in [ROOT / "WORK_MANIFEST.json", ROOT / "INPUT_MANIFEST.json", ROOT / "config.json",
                  ROOT / "Monocle2_PC.py", ROOT / "requirements.txt", ROOT / "TASK_SPEC.md"]:
            z.write(p, p.name)
        for p in sorted((ROOT / "work").glob("*.json")):
            z.write(p, p.relative_to(ROOT).as_posix())
        z.write(ROOT / "work/fit_generated.R", "work/fit_generated.R")
        z.write(target / "RESULT_MANIFEST.json", "RESULT_MANIFEST.json")
    cap = 1536 * 1024**2
    assets = []
    if archive.stat().st_size <= cap:
        assets.append(archive)
    else:
        with open(archive, "rb") as source:
            i = 1
            while True:
                first = source.read(8 * 1024**2)
                if not first:
                    break
                part = target / f"results.zip.part{i:03d}"
                with open(part, "wb") as out:
                    out.write(first)
                    written = len(first)
                    while written < cap:
                        block = source.read(min(8 * 1024**2, cap - written))
                        if not block:
                            break
                        out.write(block)
                        written += len(block)
                assets.append(part)
                i += 1
    info = {"label": label, "archive_sha256": digest(archive), "archive_bytes": archive.stat().st_size,
            "assets": [{"name": p.name, "bytes": p.stat().st_size, "sha256": digest(p)} for p in assets],
            "transfer_splitting_only": True,
            "scientific_acceptance": "Not inferred from packaging; inspect RUN_SUMMARY and acceptance records"}
    save_json(target / "RETURN_MANIFEST.json", info)
    print(json.dumps(info, indent=2), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("verify-input")
    sub.add_parser("prepare")
    f = sub.add_parser("fit")
    f.add_argument("--rscript", required=True)
    f.add_argument("--method", choices=["Monocle2", "Slingshot"], default="Monocle2")
    r = sub.add_parser("pack-return")
    r.add_argument("--label", required=True)
    args = p.parse_args()
    if args.command == "verify-input":
        verify_data()
    elif args.command == "prepare":
        prepare()
    elif args.command == "fit":
        fit(args.rscript, args.method)
    else:
        pack_return(args.label)


if __name__ == "__main__":
    main()
