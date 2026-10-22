"""Preregistered within-specimen spline model and linkage-cluster CR1 tests.

Observations are equally weighted specimen-bin expression means. These tests
are conditional on a data-derived trajectory, not independent biological proof.
"""
import numpy as np
import pandas as pd
from scipy.interpolate import BSpline
from scipy.stats import f as f_distribution


def spline_basis(t):
    t = np.asarray(t, dtype=float)
    if not np.isfinite(t).all() or t.min() < 0 or t.max() > 1:
        raise ValueError('Spline times must be finite in [0,1]')
    knots = np.array([0, 0, 0, 0, 1/3, 2/3, 1, 1, 1, 1], dtype=float)
    return BSpline.design_matrix(t, knots, 3).toarray()[:, 1:]


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    result = np.full_like(p, np.nan)
    finite = np.isfinite(p)
    if ((p[finite] < 0) | (p[finite] > 1)).any():
        raise ValueError('Invalid probability')
    indices = np.flatnonzero(finite)
    order = np.argsort(p[finite], kind='stable')
    values = p[finite][order] * len(order) / np.arange(1, len(order) + 1)
    corrected = np.minimum.accumulate(values[::-1])[::-1]
    result[indices[order]] = np.minimum(corrected, 1)
    return result


def prepare_design(observations, bins=30):
    """Encode specimen fixed intercepts, with one unit of weight per specimen."""
    required = ['specimen_id', 'cluster_unit', 'bin']
    if observations[required].isna().any().any() or observations.duplicated(['specimen_id', 'bin']).any():
        raise ValueError('Need one row for each observed specimen/bin')
    if not observations.bin.between(0, bins - 1).all():
        raise ValueError('Invalid fixed bin')
    counts = observations.groupby('specimen_id').size()
    # Singleton specimens have no within-specimen time information. They remain
    # in descriptive tables; omit algebraic zero contributions from inference.
    keep = observations.specimen_id.map(counts).ge(2).to_numpy()
    obs = observations.loc[keep].reset_index(drop=True)
    unit, units = pd.factorize(obs.specimen_id, sort=False)
    cluster, clusters = pd.factorize(obs.cluster_unit, sort=False)
    if obs.groupby('specimen_id').cluster_unit.nunique().gt(1).any():
        raise ValueError('A specimen cannot change linkage cluster between bins')
    n, u, g, k = len(obs), len(units), len(clusters), 5
    if g <= k + 1 or n <= u + k:
        raise ValueError('Insufficient linkage clusters or residual degrees of freedom')
    x = spline_basis((obs.bin.to_numpy() + .5) / bins)
    nu = np.bincount(unit, minlength=u)
    means = np.zeros((u, k))
    np.add.at(means, unit, x)
    means /= nu[:, None]
    x = x - means[unit]
    weights = 1 / nu[unit]
    information = x.T @ (weights[:, None] * x)
    if np.linalg.matrix_rank(information) != k:
        raise ValueError('Within-specimen time basis is rank deficient')
    bread = np.linalg.inv(information)
    correction = g / (g - 1) * (n - 1) / (n - u - k)
    return dict(keep=keep, obs=obs, x=x, weights=weights, unit=unit,
                units=units, unit_counts=nu, cluster=cluster, clusters=clusters,
                bread=bread, correction=correction, residual_df=n-u-k,
                bin_grid=(np.arange(bins)+.5)/bins)


def fit_dynamic_genes(observations, expression, bins=30, chunk=128):
    """Vectorize genes only; every qualifying statistical unit stays in each fit."""
    y = np.asarray(expression, dtype=float)
    if y.ndim != 2 or len(y) != len(observations) or not np.isfinite(y).all():
        raise ValueError('Expected finite specimen-bin by gene expression')
    d = prepare_design(observations, bins)
    y = y[d['keep']]
    x, w = d['x'], d['weights']
    n, genes = y.shape
    k, g, u = x.shape[1], len(d['clusters']), len(d['units'])
    beta = np.empty((k, genes))
    statistic = np.full(genes, np.nan)
    p = np.full(genes, np.nan)
    covariance_rank = np.zeros(genes, dtype=int)
    for lo in range(0, genes, chunk):
        hi = min(lo + chunk, genes)
        values = y[:, lo:hi]
        unit_means = np.zeros((u, hi-lo))
        np.add.at(unit_means, d['unit'], values)
        unit_means /= d['unit_counts'][:, None]
        centered = values - unit_means[d['unit']]
        b = d['bread'] @ (x.T @ (w[:, None] * centered))
        beta[:, lo:hi] = b
        residual = centered - x @ b
        scores = np.zeros((g, k, hi-lo))
        np.add.at(scores, d['cluster'], (w[:, None] * x)[:, :, None] * residual[:, None, :])
        transformed = np.einsum('ij,cjg->cig', d['bread'], scores)
        covariance = np.einsum('cig,cjg->gij', transformed, transformed) * d['correction']
        eigenvalues = np.linalg.eigvalsh(covariance)
        maximum = np.maximum(eigenvalues[:, -1], 0)
        rank = (eigenvalues > np.maximum(maximum[:, None] * 1e-10, 1e-20)).sum(axis=1)
        covariance_rank[lo:hi] = rank
        assessable = rank == k
        if assessable.any():
            coefficients = b[:, assessable].T
            solved = np.linalg.solve(covariance[assessable], coefficients[..., None])[..., 0]
            values_f = np.maximum(0, np.einsum('gi,gi->g', coefficients, solved) / k)
            indices = np.arange(lo, hi)[assessable]
            statistic[indices] = values_f
            p[indices] = f_distribution.sf(values_f, k, g-1)
    curves = spline_basis(d['bin_grid']) @ beta
    return dict(beta=beta, adjusted_curves=curves, statistic=statistic, p=p,
                covariance_rank=covariance_rank, observations=n,
                specimens=u, linkage_clusters=g, numerator_df=k, denominator_df=g-1,
                residual_df=d['residual_df'], CR1_correction=d['correction'],
                discarded_singleton_specimens=int(observations.specimen_id.nunique()-u),
                interpretation='Conditional within-specimen association with inferred pseudotime; not longitudinal or causal')
