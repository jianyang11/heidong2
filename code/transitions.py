"""State classification and transition-luminosity extraction per outburst.

Hardness thresholds are data-driven: the global (all-source) hardness
distribution per instrument is strongly bimodal (hard/soft branches); the
valley defines the transition hardness H_t. Transitions are located as
persistent crossings of H_t within each outburst; the transition luminosity is
interpolated at the crossing time. Errors from light-curve bootstrap.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils import (DATA_PROC, FIGS, RESULTS, TABLES, CRAB_FLUX, KPC_CM,
                   l_edd, r_s_km, get_logger)

log = get_logger("transitions")

BAND_FLUX = {"maxi": CRAB_FLUX["maxi_2_20"], "asm": CRAB_FLUX["asm_1p5_12"]}
# bolometric correction from monitor band to 0.1-100 keV, hard state ~ x3,
# at transition state we adopt 2.5 +- 1 (Migliari & Fender 06; Maccarone 03 used 7.3/3)
BOLCORR = 2.5


def hr_series(df, inst):
    g = df[(df["inst"] == inst) & np.isfinite(df["hr"]) & (df["hr"] > 0)
           & np.isfinite(df["hr_err"]) & (df["hr_err"] > 0)]
    return g[g["hr_err"] < 0.7 * g["hr"]]


def global_thresholds(sample):
    """Find valley of log-HR distribution per instrument."""
    out = {}
    for inst in ("maxi", "asm"):
        allhr = []
        for name in sample["name"]:
            p = os.path.join(DATA_PROC, f"lc_{name.replace(' ', '_')}.csv")
            if not os.path.exists(p):
                continue
            df = pd.read_csv(p)
            allhr.append(np.log10(hr_series(df, inst)["hr"]))
        v = np.concatenate(allhr)
        hist, edges = np.histogram(v, bins=60, range=(-1.5, 1.0))
        # smooth
        k = np.convolve(hist, np.ones(5) / 5, mode="same")
        c = 0.5 * (edges[1:] + edges[:-1])
        # valley: minimum between the two highest maxima
        imax1 = np.argmax(k)
        # second peak on the other side
        mask = np.abs(c - c[imax1]) > 0.35
        imax2 = np.argmax(np.where(mask, k, -1))
        lo, hi = sorted((imax1, imax2))
        ival = lo + np.argmin(k[lo:hi + 1])
        out[inst] = 10**c[ival]
        log.info("threshold %s: HR_t=%.3f (peaks at %.2f, %.2f)", inst,
                 out[inst], 10**c[imax1], 10**c[imax2])
    return out


def crossings(mjd, hr, ht, min_persist=3):
    """Return list of (idx, direction) persistent threshold crossings.
    direction -1: hard->soft (rise transition), +1: soft->hard (decay)."""
    state = np.where(hr > ht, 1, -1)
    res = []
    i = 0
    n = len(state)
    while i < n - 1:
        if state[i + 1] != state[i]:
            # require persistence: next min_persist points mostly new state
            nxt = state[i + 1:i + 1 + min_persist]
            prv = state[max(0, i - min_persist + 1):i + 1]
            if np.mean(nxt == state[i + 1]) >= 0.67 and np.mean(prv == state[i]) >= 0.67:
                res.append((i, -1 if state[i] == 1 else 1))
        i += 1
    return res


def interp_at(mjd, y, t):
    return np.interp(t, mjd, y)


def analyze_outburst(name, df, seg, ht, row, nboot=300):
    s, e = seg
    recs = []
    for inst in ("maxi", "asm"):
        h = hr_series(df, inst)
        h = h[(h["mjd"] >= s - 5) & (h["mjd"] <= e + 5)].sort_values("mjd")
        if len(h) < 15:
            continue
        # 5-point median smooth of HR for state logic
        hr_s = h["hr"].rolling(5, center=True, min_periods=2).median().values
        mjd = h["mjd"].values
        cs = crossings(mjd, hr_s, ht[inst])
        flux = df[(df["inst"] == inst)].sort_values("mjd")
        for idx, direction in cs:
            t0 = 0.5 * (mjd[idx] + mjd[idx + 1])
            fx = interp_at(flux["mjd"].values, flux["flux_crab"].values, t0)
            # bootstrap
            ts, fs = [], []
            rng = np.random.default_rng(42)
            for _ in range(nboot):
                pert = h["hr"].values + rng.normal(0, h["hr_err"].values)
                pert = pd.Series(pert).rolling(5, center=True, min_periods=2).median().values
                cb = crossings(mjd, pert, ht[inst])
                cands = [0.5 * (mjd[i] + mjd[i + 1]) for i, d in cb if d == direction]
                if cands:
                    tb = min(cands, key=lambda x: abs(x - t0))
                    if abs(tb - t0) < 30:
                        ts.append(tb)
                        fs.append(interp_at(flux["mjd"].values,
                                            flux["flux_crab"].values, tb))
            if len(ts) < 30:
                continue
            terr = np.std(ts)
            ferr = np.std(fs)
            D = row["D_kpc"] * KPC_CM
            lx = fx * BAND_FLUX[inst] * BOLCORR * 4 * np.pi * D**2
            recs.append({
                "name": name, "ob_start": s, "ob_end": e, "inst": inst,
                "type": "rise_h2s" if direction == -1 else "decay_s2h",
                "mjd_trans": t0, "mjd_err": terr,
                "flux_crab": fx, "flux_err": ferr,
                "Lx_bol": lx, "Lx_err": lx * (ferr / fx if fx > 0 else 0.3),
                "ledd_ratio": lx / l_edd(row["M_msun"]),
                "M": row["M_msun"], "M_err": row["M_err"],
                "Rs_km": r_s_km(row["M_msun"]),
            })
    return recs


def hid_figure(name, df, segs, trans, ht):
    safe = name.replace(" ", "_")
    n = max(1, len(segs))
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4.2), squeeze=False)
    for j, (s, e) in enumerate(segs):
        ax = axes[0][j]
        for inst, c in (("maxi", "tab:blue"), ("asm", "tab:green")):
            h = hr_series(df, inst)
            h = h[(h["mjd"] >= s) & (h["mjd"] <= e)]
            if len(h) < 5:
                continue
            sc = ax.scatter(h["hr"], h["flux_crab"], c=h["mjd"], s=6,
                            cmap="viridis", alpha=0.7)
            ax.axvline(ht[inst], color=c, ls=":", lw=1)
        tt = [t for t in trans if t["ob_start"] == s]
        for t in tt:
            ax.plot(ht[t["inst"]], t["flux_crab"], "r*" if "rise" in t["type"]
                    else "ms", ms=12, mec="k", zorder=5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Hardness")
        if j == 0:
            ax.set_ylabel("Flux (Crab)")
        ax.set_title(f"MJD {s:.0f}-{e:.0f}", fontsize=9)
    fig.suptitle(name)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, f"hid_{safe}.png"), dpi=130)
    plt.close(fig)


def main():
    sample = pd.read_csv(os.path.join(TABLES, "bh_sample.csv"))
    cat = pd.read_csv(os.path.join(DATA_PROC, "outburst_catalog.csv"))
    ht = global_thresholds(sample)
    allrecs = []
    for _, row in sample.iterrows():
        name = row["name"]
        safe = name.replace(" ", "_")
        p = os.path.join(DATA_PROC, f"lc_{safe}.csv")
        crow = cat[cat["name"] == name]
        if not os.path.exists(p) or crow.empty or not isinstance(
                crow.iloc[0]["outbursts"], str):
            continue
        df = pd.read_csv(p)
        segs = [tuple(map(float, x.split("-")[:1] + [x.split("-")[1]]))
                if False else
                (float(x.split("-")[0]), float(x.rsplit("-", 1)[-1]))
                for x in crow.iloc[0]["outbursts"].split(";")]
        recs = []
        for seg in segs:
            recs += analyze_outburst(name, df, seg, ht, row)
        log.info("%s: %d transitions in %d outbursts", name, len(recs), len(segs))
        hid_figure(name, df, segs, recs, ht)
        allrecs += recs
    out = pd.DataFrame(allrecs)
    out.to_csv(os.path.join(RESULTS, "transitions.csv"), index=False)
    log.info("wrote transitions.csv (%d rows)", len(out))


if __name__ == "__main__":
    main()
