"""Radio vs X-ray time lags for GX 339-4 outbursts (H2).

Radio: Corbel+13 (VizieR J/MNRAS/428/2500 table1, 8.6/9 GHz; J/MNRAS/431/L107
table2 for 2010-11 decay). X-ray: our cleaned monitor light curves
(BAT 15-50 keV as Compton proxy when available, else ASM/MAXI).
Lag convention: positive = radio lags X-ray.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils import DATA_RAW, DATA_PROC, FIGS, RESULTS, get_logger
from ccf import measure_lag

log = get_logger("radio_lags")


def read_vizier_tsv(path):
    lines = [l for l in open(path) if not l.startswith("#")]
    # header is first non-empty line
    lines = [l.rstrip("\n") for l in lines if l.strip()]
    hdr = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[2:]]  # skip units+dashes? handle below
    # drop dash separator rows
    rows = [r for r in rows if not set("".join(r)) <= set("- \t")]
    df = pd.DataFrame(rows, columns=hdr[:len(rows[0])] if rows else hdr)
    return df


def load_radio():
    a = read_vizier_tsv(os.path.join(DATA_RAW, "radio_gx339_corbel13_9715.tsv"))
    a = a.rename(columns=str.strip)
    out = []
    for _, r in a.iterrows():
        try:
            mjd = float(r["MJD"])
        except (ValueError, TypeError):
            continue
        # prefer 8.6 GHz, fall back to 4.8
        s, e = r.get("S8.6", ""), r.get("e_S8.6", "")
        ul = str(r.get("l_S8.6", "")).strip() == "<"
        if not str(s).strip() or ul:
            s, e = r.get("S4.8", ""), r.get("e_S4.8", "")
            ul = str(r.get("l_S4.8", "")).strip() == "<"
        try:
            s = float(s)
        except (ValueError, TypeError):
            continue
        if ul or not np.isfinite(s):
            continue
        try:
            e = float(e)
        except (ValueError, TypeError):
            e = 0.1 * s
        out.append({"mjd": mjd, "s_mjy": s, "e_mjy": max(e, 0.05 * s)})
    b = read_vizier_tsv(os.path.join(DATA_RAW,
                                     "radio_gx339_corbel13_2011decay.tsv"))
    for _, r in b.iterrows():
        try:
            mjd = float(r["MJD"])
            s = float(r["S9"]) if "S9" in b.columns else float(r["S9.0"])
            e = float(r.get("e_S9", r.get("e_S9.0", np.nan)))
        except (ValueError, TypeError, KeyError):
            continue
        out.append({"mjd": mjd, "s_mjy": s,
                    "e_mjy": e if np.isfinite(e) else 0.1 * s})
    df = pd.DataFrame(out).drop_duplicates("mjd").sort_values("mjd")
    log.info("radio epochs: %d (MJD %.0f-%.0f)", len(df), df.mjd.min(),
             df.mjd.max())
    return df


# hard-state windows per outburst: (label, t0, t1, phase)
# built from transitions.csv (rise: outburst start -> h2s transition;
# decay: s2h transition -> outburst end + 60)
def hard_windows():
    tr = pd.read_csv(os.path.join(RESULTS, "transitions.csv"))
    tr = tr[tr["name"] == "GX 339-4"]
    wins = []
    for (s, e), g in tr.groupby(["ob_start", "ob_end"]):
        rises = g[g["type"] == "rise_h2s"]["mjd_trans"]
        decays = g[g["type"] == "decay_s2h"]["mjd_trans"]
        if len(rises):
            wins.append((f"{int(s)}_rise", s - 30, float(rises.min()), "rise"))
        if len(decays):
            wins.append((f"{int(s)}_decay", float(decays.max()), e + 90, "decay"))
    # hard-only outbursts (no transitions): whole outburst
    cat = pd.read_csv(os.path.join(DATA_PROC, "outburst_catalog.csv"))
    row = cat[cat["name"] == "GX 339-4"].iloc[0]
    for seg in row["outbursts"].split(";"):
        s, e = float(seg.split("-")[0]), float(seg.rsplit("-", 1)[-1])
        if not ((tr["ob_start"] == s).any()):
            wins.append((f"{int(s)}_hardonly", s - 20, e + 30, "hardonly"))
    return sorted(wins, key=lambda w: w[1])


def main():
    radio = load_radio()
    lc = pd.read_csv(os.path.join(DATA_PROC, "lc_GX_339-4.csv"))
    recs = []
    for label, t0, t1, phase in hard_windows():
        r = radio[(radio.mjd >= t0) & (radio.mjd <= t1)]
        if len(r) < 6:
            log.info("%s: only %d radio pts, skip", label, len(r))
            continue
        # choose densest X-ray instrument in window (prefer BAT hard band)
        best = None
        for inst in ("bat", "asm", "maxi"):
            x = lc[(lc["inst"] == inst) & (lc.mjd >= t0 - 10) & (lc.mjd <= t1 + 10)]
            if len(x) >= 20 and (best is None or len(x) > len(best[1])):
                if best is None:
                    best = (inst, x)
                elif inst == "bat":
                    best = (inst, x)
        if best is None:
            log.info("%s: no X-ray coverage", label)
            continue
        inst, x = best
        res = measure_lag(x["mjd"].values, x["flux_crab"].values,
                          x["flux_err"].values, r["mjd"].values,
                          r["s_mjy"].values, r["e_mjy"].values,
                          lag_min=-25, lag_max=25, dlag=0.5, nmc=400)
        log.info("%s [%s, %d radio/%d xray pts]: lag=%.2f +%.2f/-%.2f d "
                 "rmax=%.2f q=%s", label, inst, len(r), len(x), res["lag"],
                 res["err_hi"], res["err_lo"], res["rmax"], res["quality"])
        recs.append({"window": label, "phase": phase, "inst": inst,
                     "n_radio": len(r), "n_xray": len(x),
                     "lag_d": res["lag"], "err_lo": res["err_lo"],
                     "err_hi": res["err_hi"], "rmax": res["rmax"],
                     "quality": res["quality"], "t0": t0, "t1": t1})
        # diagnostic figure
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
        axes[0].errorbar(x["mjd"], x["flux_crab"], x["flux_err"], fmt=".",
                         ms=3, label=f"X-ray ({inst})")
        ax2 = axes[0].twinx()
        ax2.errorbar(r["mjd"], r["s_mjy"], r["e_mjy"], fmt="rs", ms=4,
                     label="radio")
        axes[0].set_xlabel("MJD")
        axes[0].set_ylabel("X-ray flux (Crab)")
        ax2.set_ylabel("S (mJy)", color="r")
        axes[0].set_title(label)
        axes[1].plot(res["lags"], res["r"], "k-")
        axes[1].axvline(res["lag"], color="r", ls="--",
                        label=f"lag={res['lag']:.1f} d")
        axes[1].set_xlabel("lag (d, radio after X-ray)")
        axes[1].set_ylabel("ICCF r")
        axes[1].legend(fontsize=8)
        if len(res["cens"]):
            axes[2].hist(res["cens"], bins=30, color="gray")
        axes[2].set_xlabel("centroid lag (d)")
        axes[2].set_title(f"FR/RSS (q={res['quality']})")
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS, f"ccf_gx339_{label}.png"), dpi=130)
        plt.close(fig)
    out = pd.DataFrame(recs)
    out.to_csv(os.path.join(RESULTS, "tau_table.csv"), index=False)
    log.info("wrote tau_table.csv (%d windows)", len(out))


if __name__ == "__main__":
    main()
