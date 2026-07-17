"""Hierarchical decomposition + spin/ISCO model comparison (event-level frame).

(1) Variance decomposition: for event-level responses y (hysteresis amplitude,
    decay transition log L/L_Edd), fit y_ij = mu + u_i + e_ij with
    u_i ~ N(0, sig_s^2) (source level), e_ij ~ N(0, sig_e^2) (event level),
    marginalizing u_i analytically. Fraction f_s = sig_s^2/(sig_s^2+sig_e^2)
    quantifies how much is carried by immutable source properties (spin, M)
    vs per-event state (magnetic flux, mdot history).
(2) Jet radio-loudness xi (per-source normalization of the hard-state Lr-Lx
    track) regressed on competing spin parametrizations:
      M0: constant;  M1: xi = c + b a^2;  M2: xi = c + b Omega_H^2 (BZ);
      M3: xi = c + b log r_ISCO;  each +- phi_max covariate.
    Compared with WAIC (small-sample caveat reported).
(3) Event-level demonstration at fixed spin: GX 339-4 radio peak vs X-ray peak.
"""
import os
from itertools import product

import numpy as np
import pandas as pd
import emcee
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils import RESULTS, FIGS, get_logger

log = get_logger("hier")
OUT = []


def say(m):
    log.info(m)
    OUT.append(m)


def fit_varcomp(df, col, nstep=4000):
    groups = [g[col].values for _, g in df.groupby("name") if False] or \
             [g.values for _, g in df.groupby("name")[col]]
    groups = [g[np.isfinite(g)] for g in groups]
    groups = [g for g in groups if len(g) > 0]

    def lnpost(th):
        mu, lss, lse = th
        if not (-6 < mu < 6 and -6 < lss < 2 and -6 < lse < 2):
            return -np.inf
        ss2, se2 = np.exp(2 * lss), np.exp(2 * lse)
        lp = 0.0
        for y in groups:
            n = len(y)
            # marginal: cov = se2 I + ss2 J
            ybar = y.mean()
            dev = y - ybar
            det = se2**(n - 1) * (se2 + n * ss2)
            quad = (dev**2).sum() / se2 + n * (ybar - mu)**2 / (se2 + n * ss2)
            lp += -0.5 * (np.log(det) + quad + n * np.log(2 * np.pi))
        return lp
    rng = np.random.default_rng(11)
    p0 = np.array([np.nanmean(np.concatenate(groups)), -1.0, -1.0])
    pos = p0 + 0.05 * rng.standard_normal((24, 3))
    sam = emcee.EnsembleSampler(24, 3, lnpost)
    sam.run_mcmc(pos, nstep, progress=False)
    ch = sam.get_chain(discard=nstep // 2, flat=True)
    ss2, se2 = np.exp(2 * ch[:, 1]), np.exp(2 * ch[:, 2])
    fs = ss2 / (ss2 + se2)
    q = np.percentile(fs, [16, 50, 84])
    say(f"[varcomp {col}] n_src={len(groups)} n_evt={sum(len(g) for g in groups)} "
        f"sig_source={np.median(np.sqrt(ss2)):.3f} dex, "
        f"sig_event={np.median(np.sqrt(se2)):.3f} dex, "
        f"f_source={q[1]:.2f} (+{q[2]-q[1]:.2f}/-{q[1]-q[0]:.2f})")
    return fs


def waic_regression(x, y, ye, nstep=3000, seed=5):
    """y = c + b x, returns WAIC and posterior of b."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ye = np.asarray(ye, float)

    def lnpost(th):
        c, b, lns = th
        if not (-50 < c < 60 and -30 < b < 30 and -6 < lns < 2):
            return -np.inf
        s2 = np.exp(2 * lns) + ye**2
        return -0.5 * np.sum((y - c - b * x)**2 / s2 + np.log(2 * np.pi * s2))
    rng = np.random.default_rng(seed)
    p0 = np.array([np.mean(y), 0.0, np.log(np.std(y) + 1e-2)])
    pos = p0 + 0.05 * rng.standard_normal((24, 3))
    sam = emcee.EnsembleSampler(24, 3, lnpost)
    sam.run_mcmc(pos, nstep, progress=False)
    ch = sam.get_chain(discard=nstep // 2, flat=True)[::10]
    # pointwise log-lik
    ll = np.empty((len(ch), len(y)))
    for i, (c, b, lns) in enumerate(ch):
        s2 = np.exp(2 * lns) + ye**2
        ll[i] = -0.5 * ((y - c - b * x)**2 / s2 + np.log(2 * np.pi * s2))
    lppd = np.sum(np.log(np.mean(np.exp(ll), axis=0)))
    pwaic = np.sum(np.var(ll, axis=0))
    return -2 * (lppd - pwaic), ch[:, 1]


def main():
    ev = pd.read_csv(os.path.join(RESULTS, "events.csv"))
    src = pd.read_csv(os.path.join(RESULTS, "source_table.csv"))

    say("== (1) source-level vs event-level variance decomposition ==")
    h = ev[np.isfinite(ev["hyst_dex"])]
    fit_varcomp(h, "hyst_dex")
    d = ev[np.isfinite(ev["ledd_decay"])].copy()
    d["log_ledd_decay"] = np.log10(d["ledd_decay"])
    fit_varcomp(d, "log_ledd_decay")
    r = ev[np.isfinite(ev["ledd_rise"])].copy()
    r["log_ledd_rise"] = np.log10(r["ledd_rise"])
    fit_varcomp(r, "log_ledd_rise")

    say("== (2) jet radio-loudness xi vs spin parametrizations ==")
    s = src[np.isfinite(src["xi"])].copy()
    s["xi_err"] = np.clip(s["xi_err"], 0.15, None)  # floor for n=1 sources
    sa = s[np.isfinite(s["a_best"])]
    say(f"sources with xi: {len(s)}; with spin: {len(sa)} "
        f"({(sa['a_method']=='CF').sum()} CF, {(sa['a_method']=='refl').sum()} refl)")
    preds = {"M0_const": np.zeros(len(sa)),
             "M1_a2": sa["a_best"].values**2,
             "M2_OmegaH2": sa["omega_h"].values**2,
             "M3_logRisco": np.log10(sa["r_isco_rg"].values)}
    for lab, x in preds.items():
        w, b = waic_regression(x, sa["xi"].values, sa["xi_err"].values)
        qb = np.percentile(b, [16, 50, 84])
        say(f"  {lab}: WAIC={w:.1f} slope={qb[1]:+.2f} (+{qb[2]-qb[1]:.2f}"
            f"/-{qb[1]-qb[0]:.2f})")
    # phi covariate check
    sp = sa[np.isfinite(sa["phi_max_mid"])]
    if len(sp) >= 5:
        w, b = waic_regression(np.log10(sp["phi_max_mid"].values),
                               sp["xi"].values, sp["xi_err"].values)
        qb = np.percentile(b, [16, 50, 84])
        say(f"  M4_logphi (n={len(sp)}): WAIC={w:.1f} slope={qb[1]:+.2f} "
            f"(+{qb[2]-qb[1]:.2f}/-{qb[1]-qb[0]:.2f})")
    # spin-method sensitivity
    for meth, col in (("CF", "a_cf"), ("refl", "a_refl")):
        ss = s[np.isfinite(s[col])]
        if len(ss) >= 4:
            w, b = waic_regression(ss[col].values**2, ss["xi"].values,
                                   ss["xi_err"].values)
            qb = np.percentile(b, [16, 50, 84])
            say(f"  sensitivity a^2 ({meth} only, n={len(ss)}): "
                f"slope={qb[1]:+.2f} (+{qb[2]-qb[1]:.2f}/-{qb[1]-qb[0]:.2f})")

    say("== (3) event-level jet output at fixed spin: GX 339-4 ==")
    g = ev[(ev["name"] == "GX 339-4") & np.isfinite(ev["radio_peak_mjy"])]
    if len(g) >= 4:
        rho = np.corrcoef(np.log10(g["peak_ledd"]),
                          np.log10(g["radio_peak_mjy"]))[0, 1]
        rng_dex = np.log10(g["radio_peak_mjy"].max() / g["radio_peak_mjy"].min())
        say(f"GX 339-4: {len(g)} events, radio-peak spread {rng_dex:.2f} dex at "
            f"fixed spin; corr(log peak_ledd, log S_peak)={rho:.2f} -> "
            f"event-state drives >~{rng_dex:.1f} dex of jet variance within one source")

    # figures
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    ax = axes[0]
    for meth, c in (("CF", "tab:blue"), ("refl", "tab:red")):
        m = sa["a_method"] == meth
        ax.errorbar(sa[m]["a_best"], sa[m]["xi"], yerr=sa[m]["xi_err"],
                    fmt="o", color=c, ms=7, mec="k", capsize=2, label=meth)
    for _, rr in sa.iterrows():
        ax.annotate(rr["name"], (rr["a_best"], rr["xi"]), fontsize=6,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("spin a")
    ax.set_ylabel(r"$\xi$ = log $L_R$ at $L_X=10^{36}$ (radio loudness)")
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.errorbar(sa["r_isco_rg"], sa["xi"], yerr=sa["xi_err"], fmt="s",
                color="tab:purple", ms=7, mec="k", capsize=2)
    ax.set_xlabel(r"$r_{\rm ISCO}$ ($r_g$)")
    ax.set_ylabel(r"$\xi$")
    ax = axes[2]
    if len(g) >= 4:
        ax.plot(g["peak_ledd"], g["radio_peak_mjy"], "o", ms=8, mec="k",
                color="tab:green")
        for _, rr in g.iterrows():
            ax.annotate(f"MJD {int(rr['ob_start'])}",
                        (rr["peak_ledd"], rr["radio_peak_mjy"]), fontsize=6,
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"outburst peak $L/L_{\rm Edd}$")
        ax.set_ylabel("radio peak (mJy)")
        ax.set_title("GX 339-4: event scatter at fixed spin")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "xi_spin_isco_events.png"), dpi=140)
    plt.close(fig)

    with open(os.path.join(RESULTS, "hier_summary.txt"), "w") as f:
        f.write("\n".join(OUT) + "\n")
    log.info("wrote hier_summary.txt")


if __name__ == "__main__":
    main()
