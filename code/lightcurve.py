"""Clean monitoring light curves, build daily products, HIDs, outburst catalog.

Products per source (data/processed/):
  lc_<src>.csv : MJD, inst, flux_crab (band-summed), hr (hard/soft), errors
Figures (output/figures/): lc_<src>.png overview.
"""
import glob
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

from utils import (DATA_RAW, DATA_PROC, FIGS, TABLES, CRAB_MAXI, CRAB_BAT,
                   CRAB_ASM, get_logger)

log = get_logger("lightcurve")

MJDREF_ASM = 49353.000696574074  # RXTE MJDREF


def load_maxi(path):
    # cols: MJD F(2-20) e F(2-4) e F(4-10) e F(10-20) e
    d = np.loadtxt(path)
    df = pd.DataFrame(d, columns=["mjd", "f220", "e220", "f24", "e24",
                                  "f410", "e410", "f1020", "e1020"])
    df["flux_crab"] = df["f220"] / CRAB_MAXI["2-20"]
    df["flux_err"] = df["e220"] / CRAB_MAXI["2-20"]
    # hardness: 4-10 / 2-4 in Crab units of each band
    with np.errstate(divide="ignore", invalid="ignore"):
        df["hr"] = df["f410"] / df["f24"]
        df["hr_err"] = df["hr"] * np.sqrt((df["e410"] / df["f410"])**2 +
                                          (df["e24"] / df["f24"])**2)
    df["inst"] = "maxi"
    return df[["mjd", "inst", "flux_crab", "flux_err", "hr", "hr_err"]]


def load_bat(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            try:
                rows.append((float(p[0]), float(p[1]), float(p[2])))
            except (ValueError, IndexError):
                continue
    df = pd.DataFrame(rows, columns=["mjd", "rate", "err"])
    df["flux_crab"] = df["rate"] / CRAB_BAT
    df["flux_err"] = df["err"] / CRAB_BAT
    df["hr"] = np.nan
    df["hr_err"] = np.nan
    df["inst"] = "bat"
    return df[["mjd", "inst", "flux_crab", "flux_err", "hr", "hr_err"]]


def _read_asm_fits(path):
    with fits.open(path) as h:
        d = h[1].data
        mjd = np.asarray(d["TIME"], dtype="=f8") + 49353.0  # TIMEUNIT is days
        band = np.char.strip(d["BAND"].astype(str))
        rate = np.asarray(d["RATE"], dtype="=f8")
        err = np.asarray(d["ERROR"], dtype="=f8")
    return mjd, band, rate, err


def load_asm(path):
    mjd, band, rate, err = _read_asm_fits(path)
    rows = {"SUM": pd.DataFrame({"mjd": mjd, "rate": rate, "err": err})}
    colpath = path.replace("asm_", "asmcol_").replace(".lc", ".col")
    if os.path.exists(colpath):
        cm, cb, cr, ce = _read_asm_fits(colpath)
        for b in ("A", "B", "C"):
            m = cb == b
            rows[b] = pd.DataFrame({"mjd": cm[m], "rate": cr[m], "err": ce[m]})
    else:
        for b in ("A", "B", "C"):
            rows[b] = pd.DataFrame({"mjd": [], "rate": [], "err": []})
    # daily bins
    def daily(df):
        df = df[np.isfinite(df["rate"]) & np.isfinite(df["err"]) & (df["err"] > 0)].copy()
        df["day"] = np.floor(df["mjd"])
        df["w"] = 1 / df["err"]**2
        df["rw"] = df["rate"] * df["w"]
        g = df.groupby("day")[["rw", "w"]].sum()
        out = pd.DataFrame({"mjd": g.index.values + 0.5,
                            "rate": g["rw"] / g["w"],
                            "err": 1 / np.sqrt(g["w"])}).reset_index(drop=True)
        return out
    s, a, c = daily(rows["SUM"]), daily(rows["A"]), daily(rows["C"])
    df = s.merge(a, on="mjd", suffixes=("", "_a")).merge(
        c, on="mjd", suffixes=("", "_c"))
    df["flux_crab"] = df["rate"] / CRAB_ASM["sum"]
    df["flux_err"] = df["err"] / CRAB_ASM["sum"]
    with np.errstate(divide="ignore", invalid="ignore"):
        df["hr"] = df["rate_c"] / df["rate_a"]
        df["hr_err"] = np.abs(df["hr"]) * np.sqrt(
            (df["err_c"] / df["rate_c"])**2 + (df["err_a"] / df["rate_a"])**2)
    df["inst"] = "asm"
    return df[["mjd", "inst", "flux_crab", "flux_err", "hr", "hr_err"]]


def clean(df, snr=2.0):
    df = df[np.isfinite(df["flux_crab"]) & np.isfinite(df["flux_err"])]
    df = df[df["flux_err"] > 0]
    df = df[df["flux_crab"] / df["flux_err"] > snr]
    df = df[df["flux_crab"] > 0]
    return df.sort_values("mjd").reset_index(drop=True)


# per-source overrides for faint/hard-only or short outbursts
OB_OVERRIDES = {
    "XTE J1118+480": {"min_peak": 0.025, "thresh_crab": 0.02},
    "V404 Cyg": {"min_days": 8, "insts": ("maxi", "asm", "bat")},
    "MAXI J1659-152": {"min_days": 10},
}


def find_outbursts(df, thresh_crab=0.03, min_days=15, gap_days=30,
                   min_peak=0.1, min_pts=10, insts=("maxi", "asm")):
    """Segments where flux > thresh (Crab) with persistence; merge close segments."""
    d = df[df["inst"].isin(insts)]
    if len(d) < 30:
        d = df
    active = d[d["flux_crab"] > thresh_crab]
    if len(active) == 0:
        return []
    mjds = active["mjd"].values
    segs = []
    start = prev = mjds[0]
    for m in mjds[1:]:
        if m - prev > gap_days:
            segs.append((start, prev))
            start = m
        prev = m
    segs.append((start, prev))
    out = []
    for s, e in segs:
        if e - s < min_days:
            continue
        seg = active[(active["mjd"] >= s) & (active["mjd"] <= e)]
        if len(seg) >= min_pts and seg["flux_crab"].max() >= min_peak:
            out.append((s, e))
    return out


def process_source(name):
    safe = name.replace(" ", "_")
    parts = []
    for pat, loader in ((f"maxi_{safe}.dat", load_maxi),
                        (f"bat_{safe}.dat", load_bat),
                        (f"asm_{safe}.lc", load_asm)):
        p = os.path.join(DATA_RAW, pat)
        if os.path.exists(p):
            try:
                parts.append(clean(loader(p)))
                log.info("%s: loaded %s (%d pts)", name, pat, len(parts[-1]))
            except Exception as e:
                log.error("%s: %s failed: %s", name, pat, e)
    if not parts:
        log.warning("%s: no data", name)
        return None
    df = pd.concat(parts).sort_values("mjd").reset_index(drop=True)
    df.to_csv(os.path.join(DATA_PROC, f"lc_{safe}.csv"), index=False)
    obs = find_outbursts(df, **OB_OVERRIDES.get(name, {}))
    # overview figure
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    colors = {"maxi": "tab:blue", "asm": "tab:green", "bat": "tab:red"}
    for inst, g in df.groupby("inst"):
        axes[0].errorbar(g["mjd"], g["flux_crab"], g["flux_err"], fmt=".",
                         ms=2, alpha=0.5, color=colors[inst], label=inst.upper())
        gh = g[np.isfinite(g["hr"]) & np.isfinite(g["hr_err"]) & (g["hr"] > 0)
               & (g["hr_err"] > 0) & (g["hr_err"] < 0.5 * np.abs(g["hr"]))]
        axes[1].errorbar(gh["mjd"], gh["hr"], gh["hr_err"], fmt=".", ms=2,
                         alpha=0.5, color=colors[inst])
    for s, e in obs:
        axes[0].axvspan(s, e, color="gold", alpha=0.2)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Flux (Crab)")
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].set_title(f"{name}  ({len(obs)} outburst segments)")
    hr_ok = df[np.isfinite(df["hr"]) & (df["hr"] > 0)]
    if len(hr_ok) > 5:
        axes[1].set_yscale("log")
    axes[1].set_ylabel("Hardness")
    axes[1].set_xlabel("MJD")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, f"lc_{safe}.png"), dpi=130)
    plt.close(fig)
    return {"name": name, "npts": len(df), "outbursts": obs}


def main():
    sample = pd.read_csv(os.path.join(TABLES, "bh_sample.csv"))
    inv = []
    for name in sample["name"]:
        r = process_source(name)
        if r:
            inv.append({"name": name, "npts": r["npts"],
                        "n_outbursts": len(r["outbursts"]),
                        "outbursts": ";".join(f"{s:.0f}-{e:.0f}" for s, e in r["outbursts"])})
    pd.DataFrame(inv).to_csv(os.path.join(DATA_PROC, "outburst_catalog.csv"), index=False)
    log.info("wrote outburst_catalog.csv (%d sources)", len(inv))


if __name__ == "__main__":
    main()
