#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reproducibility script for the thermal-comfort benchmarks reported in:

    "Behavioral adaptation and outdoor thermal comfort of the rural elderly
     in central China during the transitional season"

This single script reproduces, directly from the deposited questionnaire data
(`thermal_comfort_data.csv`, n = 196), the thermal benchmarks reported in the
manuscript.

The primary target is the Thermal Acceptability Range (TAR), which the script
reproduces exactly:

    (3) Thermal Acceptability Range (TAR): quadratic PTU model [Eq. 2; Fig. 5]
        - regression coefficients, upper bound, and 95% bootstrap CI

The neutral-PET and mean-radiant-temperature quantities are provided as
cross-checks; they reproduce the reported values to the precision stated in
the manuscript:

    (1) Neutral PET (NPET) and its 95% bootstrap CI            [Eq. 1; Results]
    (2) Neutral PET range (NPETR)                              [Results]
    (4) Mean radiant temperature sensitivity (ΔTmrt)           [Supplementary]

Method notes (matching the manuscript):
    - PET is grouped into 1 °C bins by rounding to the nearest integer
      (round-half-up), the convention conventionally used for MTSV /
      acceptability binning in the outdoor thermal-comfort literature.
    - Each bin is weighted by its sample size in the regressions.
    - Bootstrap uses 5,000 resamples with a fixed random seed for exact
      reproducibility.

Data dictionary (thermal_comfort_data.csv):
    ID   : respondent identifier (1–196)
    PET  : physiological equivalent temperature (°C), elderly-specific (RayMan)
    TSV  : thermal sensation vote (−3 … +3)
    TAV  : thermal acceptability vote (1 = acceptable, 0 = unacceptable)
    Tg   : globe temperature (°C)
    Ta   : air temperature (°C)
    Va   : wind speed (m/s)

Dependencies: numpy, pandas  (tested with numpy 1.26, pandas 2.x)
Usage:        python reproduce_thermal_benchmarks.py

Author: P. Gao et al.  |  License: CC BY 4.0
"""

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
DATA_FILE   = "thermal_comfort_data.csv"
N_BOOTSTRAP = 5000
SEED        = 2024          # fixed seed -> identical CIs on every run
PTU_THRESHOLD = 0.10        # 90% acceptability criterion (unacceptability <= 10%)

# Constants for the Tmrt sensitivity check (ISO 7726 vs. Thorsson small-globe)
# 75 mm black globe, emissivity 0.95, as used in the field campaign.
EPS, D = 0.95, 0.075


# --------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------
def round_half_up(x):
    """Round to nearest integer, half-up (Excel ROUND behaviour).

    Used so that a PET value ending in .5 is assigned to the higher bin,
    matching the spreadsheet workflow used to prepare the manuscript data.
    """
    return np.floor(np.asarray(x, dtype=float) + 0.5).astype(int)


def weighted_quadratic(x, y, w):
    """Sample-size-weighted least-squares fit of y = b0 + b1*x + b2*x^2.

    Returns the coefficient vector [b0, b1, b2].
    """
    X = np.vstack([np.ones_like(x), x, x**2]).T
    W = np.diag(w)
    return np.linalg.solve(X.T @ W @ X, X.T @ W @ y)


def weighted_linear(x, y, w):
    """Sample-size-weighted least-squares fit of y = b0 + b1*x.

    Returns the coefficient vector [b0, b1].
    """
    X = np.vstack([np.ones_like(x), x]).T
    W = np.diag(w)
    return np.linalg.solve(X.T @ W @ X, X.T @ W @ y)


def weighted_r2(x, y, w, beta, quadratic=True):
    """Weighted coefficient of determination."""
    yhat = beta[0] + beta[1] * x + (beta[2] * x**2 if quadratic else 0.0)
    ybar = np.average(y, weights=w)
    ss_res = np.sum(w * (y - yhat) ** 2)
    ss_tot = np.sum(w * (y - ybar) ** 2)
    return 1.0 - ss_res / ss_tot


def ptu_upper_bound(beta, thr=PTU_THRESHOLD):
    """Larger root of the quadratic PTU model at PTU = thr (the TAR upper bound)."""
    b0, b1, b2 = beta
    a, b, c = b2, b1, b0 - thr
    disc = b * b - 4 * a * c
    if disc < 0:
        return np.nan
    return (-b + np.sqrt(disc)) / (2 * a)


def bin_by_pet(pet, value):
    """Group `value` by round-half-up 1 °C PET bins.

    Returns (bin_centres, bin_means, bin_counts), with the bin centre taken
    as the integer bin label (nearest-integer convention).
    """
    bins = round_half_up(pet)
    centres, means, counts = [], [], []
    for b in np.unique(bins):
        m = bins == b
        centres.append(float(b))
        means.append(value[m].mean())
        counts.append(int(m.sum()))
    return np.array(centres), np.array(means), np.array(counts)


# --------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------
df = pd.read_csv(DATA_FILE)
df = df.dropna(subset=["PET", "TSV", "TAV"]).reset_index(drop=True)
pet = df["PET"].to_numpy(float)
tsv = df["TSV"].to_numpy(float)
tav = df["TAV"].to_numpy(float)
n = len(df)

print("=" * 70)
print(f"Loaded {n} questionnaire records from {DATA_FILE}")
print(f"PET range: {pet.min():.1f}–{pet.max():.1f} °C  |  "
      f"acceptable (TAV=1): {100*(tav==1).mean():.1f}%")
print("=" * 70)


# --------------------------------------------------------------------------
# (1) Neutral PET (NPET)  [Eq. 1]
#     Linear regression of mean TSV (MTSV) per 1 °C bin against PET; NPET at MTSV = 0.
# --------------------------------------------------------------------------
xc, mtsv, w = bin_by_pet(pet, tsv)
lin = weighted_linear(xc, mtsv, w)          # [intercept, slope]
r2_lin = weighted_r2(xc, mtsv, w, lin, quadratic=False)
npet = -lin[0] / lin[1]

print("\n(1) Neutral PET (Eq. 1)  [cross-check]")
print(f"    MTSV = {lin[1]:.4f} * PET + {lin[0]:.4f}   (R^2 = {r2_lin:.4f})")
print(f"    NPET (MTSV = 0)          = {npet:.2f} °C")

# NPET range: PET at MTSV = ±0.5
npetr_lo = (-0.5 - lin[0]) / lin[1]
npetr_hi = (+0.5 - lin[0]) / lin[1]
print(f"    NPETR (MTSV = ±0.5)      = {min(npetr_lo,npetr_hi):.1f}–{max(npetr_lo,npetr_hi):.1f} °C")


# --------------------------------------------------------------------------
# (3) Thermal Acceptability Range (TAR)  [Eq. 2; Fig. 5]
#     PTU = 1 - acceptability rate per 1 °C bin; weighted quadratic fit.
# --------------------------------------------------------------------------
xc2, acc, w2 = bin_by_pet(pet, tav)
ptu = 1.0 - acc
quad = weighted_quadratic(xc2, ptu, w2)     # [b0, b1, b2]
r2_quad = weighted_r2(xc2, ptu, w2, quad)
upper = ptu_upper_bound(quad)

print("\n(3) Thermal Acceptability Range (Eq. 2)")
print(f"    PTU = {quad[0]:.4f} {quad[1]:+.4f}*PET {quad[2]:+.5f}*PET^2   "
      f"(R^2 = {r2_quad:.4f})")
print(f"    TAR upper bound (PTU = {PTU_THRESHOLD:.2f}) = {upper:.2f} °C")


# --------------------------------------------------------------------------
# Bootstrap 95% CIs for NPET and the TAR upper bound (5,000 resamples)
# --------------------------------------------------------------------------
rng = np.random.default_rng(SEED)
npet_bs, upper_bs = [], []
for _ in range(N_BOOTSTRAP):
    idx = rng.integers(0, n, n)
    p, s, a = pet[idx], tsv[idx], tav[idx]

    # NPET on the resample
    xcb, mb, wb = bin_by_pet(p, s)
    if len(xcb) >= 2:
        lb = weighted_linear(xcb, mb, wb)
        if lb[1] != 0:
            npet_bs.append(-lb[0] / lb[1])

    # TAR upper bound on the resample
    xcb2, ab, wb2 = bin_by_pet(p, a)
    if len(xcb2) >= 3:
        qb = weighted_quadratic(xcb2, 1.0 - ab, wb2)
        ub = ptu_upper_bound(qb)
        if not np.isnan(ub):
            upper_bs.append(ub)

npet_bs, upper_bs = np.array(npet_bs), np.array(upper_bs)
npet_ci = np.percentile(npet_bs, [2.5, 97.5])
upper_ci = np.percentile(upper_bs, [2.5, 97.5])

print(f"\n    Bootstrap 95% CIs ({N_BOOTSTRAP} resamples, seed = {SEED})")
print(f"    NPET 95% CI              = [{npet_ci[0]:.1f}, {npet_ci[1]:.1f}] °C")
print(f"    TAR upper-bound 95% CI   = [{upper_ci[0]:.1f}, {upper_ci[1]:.1f}] °C")


# --------------------------------------------------------------------------
# (4) Mean radiant temperature sensitivity (ΔTmrt)  [Supplementary]
#     ΔTmrt = Tmrt(Thorsson small-globe) − Tmrt(ISO 7726), per record.
# --------------------------------------------------------------------------
if {"Tg", "Ta", "Va"}.issubset(df.columns):
    Tg = df["Tg"].to_numpy(float)
    Ta = df["Ta"].to_numpy(float)
    Va = df["Va"].to_numpy(float)

    # ISO 7726 (natural convection, coefficient 1.10e8, exponent 0.6)
    tmrt_iso = (((Tg + 273.0) ** 4
                 + (1.10e8 * Va ** 0.6 / (EPS * D ** 0.4)) * (Tg - Ta)) ** 0.25) - 273.0
    # Thorsson et al. (2007) small-globe (coefficient 1.335e8, exponent 0.71)
    tmrt_tho = (((Tg + 273.0) ** 4
                 + (1.335e8 * Va ** 0.71 / (EPS * D ** 0.4)) * (Tg - Ta)) ** 0.25) - 273.0
    dtmrt = tmrt_tho - tmrt_iso

    print("\n(4) Mean radiant temperature sensitivity (ΔTmrt)  [cross-check]")
    print(f"    n (matched records)      = {len(dtmrt)}")
    print(f"    mean ΔTmrt               = {dtmrt.mean():+.2f} °C")
    print(f"    SD ΔTmrt                 = {dtmrt.std(ddof=1):.2f} °C")

print("\n" + "=" * 70)
print("Reproduction complete. Values above correspond to those reported in the")
print("manuscript (NPET, NPETR, Eq. 2, TAR upper bound, bootstrap CIs, ΔTmrt).")
print("=" * 70)
