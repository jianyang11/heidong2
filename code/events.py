"""Build the unified event-level sample (one row per outburst event).

Merges: outburst catalog (peak L/L_Edd from cleaned light curves), transition
luminosities (rise/decay), hysteresis amplitude, per-event GX 339-4 radio
peaks (Corbel+13 measurements), source-level spin table (CF & reflection),
r_ISCO, Omega_H, BZ efficiency factor, and per-source jet-normalization xi
from the hard-state Lr-Lx database.

xi := median residual of log Lr about the global hard-state correlation
log Lr = xi_global + 0.6 (log Lx - 36); a distance-independent* per-source
jet radio-loudness at fixed X-ray luminosity (*distance enters both L's;
documented as caveat).
"""
import os

import numpy as np
import pandas as pd

from utils import (BASE, DATA_PROC, RESULTS, TABLES, CRAB_FLUX, KPC_CM,
                   l_edd, get_logger)
from transitions import BAND_FLUX, BOLCORR
from phi_jet import ALIASES
from radio_lags import load_radio

log = get_logger("events")


def r_isco_rg(a):
    """Bardeen, Press & Teukolsky 1972; prograde. In units of r_g."""
    z1 = 1 + (1 - a**2)**(1 / 3) * ((1 + a)**(1 / 3) + (1 - a)**(1 / 3))
    z2 = np.sqrt(3 * a**2 + z1**2)
    return 3 + z2 - np.sqrt((3 - z1) * (3 + z1 + 2 * z2))


def omega_h(a):
    """Dimensionless horizon angular frequency a/(2 r_H/r_g)."""
    return a / (2 * (1 + np.sqrt(1 - min(a, 0.998)**2)))


def jet_norm():
    db = pd.read_csv(os.path.join(BASE, "data", "XRB-LrLx_pub", "data",
                                  "lrlx_data_BHs.csv"))
    db["name_std"] = db["Name"].map(lambda n: ALIASES.get(str(n).strip()))
    db = db[db["name_std"].notna() & db["uplim"].isna()
            & (db.Lr > 0) & (db.Lx > 0)]
    db = db[(db.Lx > 1e33) & (db.Lx < 3e37)]  # hard-state track, no quiescence
    res = np.log10(db["Lr"]) - 0.6 * (np.log10(db["Lx"]) - 36)
    out = []
    rng = np.random.default_rng(3)
    for name, g in res.groupby(db["name_std"]):
        boots = [np.median(rng.choice(g, len(g))) for _ in range(500)]
        out.append({"name": name, "xi": np.median(g), "xi_err": np.std(boots),
                    "n_lrlx": len(g)})
    return pd.DataFrame(out)


def main():
    sample = pd.read_csv(os.path.join(TABLES, "bh_sample.csv"))
    spin = pd.read_csv(os.path.join(TABLES, "spin_table.csv"))
    cat = pd.read_csv(os.path.join(DATA_PROC, "outburst_catalog.csv"))
    tr = pd.read_csv(os.path.join(RESULTS, "transitions.csv"))
    veto = pd.read_csv(os.path.join(TABLES, "veto_windows.csv"))
    radio = load_radio()

    rows = []
    for _, s in cat.iterrows():
        name = s["name"]
        srow = sample[sample["name"] == name].iloc[0]
        lc = pd.read_csv(os.path.join(
            DATA_PROC, f"lc_{name.replace(' ', '_')}.csv"))
        for seg in str(s["outbursts"]).split(";"):
            p = seg.strip().rsplit("-", 1)
            if len(p) != 2 or not p[0]:
                continue
            t0, t1 = float(p[0]), float(p[1])
            # veto drops
            v = veto[(veto["name"] == name) & (veto["action"] == "drop")
                     & (veto["mjd_lo"] <= t0) & (veto["mjd_hi"] >= t1)]
            if len(v):
                continue
            m = (lc.mjd >= t0) & (lc.mjd <= t1) & lc["inst"].isin(["maxi", "asm"])
            if not m.any():
                continue
            sub = lc[m]
            irow = sub["flux_crab"].idxmax()
            pk = sub.loc[irow, "flux_crab"]
            inst_pk = sub.loc[irow, "inst"]
            # same conversion as transitions.py
            lx_pk = (pk * BAND_FLUX[inst_pk] * BOLCORR * 4 * np.pi
                     * (srow["D_kpc"] * KPC_CM)**2)
            ledd_pk = lx_pk / l_edd(srow["M_msun"])
            g = tr[(tr["name"] == name) & (tr["ob_start"] == t0)]
            r = g[g["type"] == "rise_h2s"]["ledd_ratio"]
            d = g[g["type"] == "decay_s2h"]["ledd_ratio"]
            row = {"name": name, "ob_start": t0, "ob_end": t1,
                   "peak_crab": pk, "peak_ledd": ledd_pk,
                   "ledd_rise": r.max() if len(r) else np.nan,
                   "ledd_decay": d.min() if len(d) else np.nan,
                   "hyst_dex": (np.log10(r.max() / d.min())
                                if len(r) and len(d) else np.nan),
                   "radio_peak_mjy": np.nan, "radio_ref": ""}
            if name == "GX 339-4":
                rm = radio[(radio.mjd >= t0 - 30) & (radio.mjd <= t1 + 90)]
                if len(rm) >= 3:
                    row["radio_peak_mjy"] = rm["s_mjy"].max()
                    row["radio_ref"] = "Corbel+13 (VizieR); >=3 epochs in event"
            rows.append(row)
    ev = pd.DataFrame(rows)
    # attach source-level
    src = sample.merge(spin, on="name", how="left").merge(jet_norm(), on="name",
                                                          how="left")
    phi = pd.read_csv(os.path.join(RESULTS, "phi_jet_per_source.csv"))
    src = src.merge(phi[["name", "phi_max_mid", "phi_max_low", "phi_max_high"]],
                    on="name", how="left")
    a_best = src["a_cf"].where(src["a_cf"].notna(), src["a_refl"])
    src["a_best"] = a_best
    src["a_method"] = np.where(src["a_cf"].notna(), "CF",
                               np.where(src["a_refl"].notna(), "refl", "none"))
    src["r_isco_rg"] = [r_isco_rg(a) if np.isfinite(a) else np.nan
                        for a in a_best]
    src["omega_h"] = [omega_h(a) if np.isfinite(a) else np.nan for a in a_best]
    src["eta_bz"] = src["omega_h"]**2
    keep = ["name", "M_msun", "M_err", "D_kpc", "D_err", "Porb_d", "tier",
            "a_cf", "a_cf_err", "a_refl", "a_refl_err", "a_best", "a_method",
            "r_isco_rg", "omega_h", "eta_bz", "xi", "xi_err", "n_lrlx",
            "phi_max_mid", "phi_max_low", "phi_max_high"]
    src[keep].to_csv(os.path.join(RESULTS, "source_table.csv"), index=False)
    ev = ev.merge(src[keep], on="name", how="left")
    ev.to_csv(os.path.join(RESULTS, "events.csv"), index=False)
    log.info("events: %d rows, %d sources; with hysteresis: %d; with radio "
             "peak: %d", len(ev), ev["name"].nunique(),
             ev["hyst_dex"].notna().sum(), ev["radio_peak_mjy"].notna().sum())
    log.info("source xi:\n%s", src[["name", "xi", "xi_err", "n_lrlx",
                                    "a_best", "a_method"]].dropna(
        subset=["xi"]).to_string())


if __name__ == "__main__":
    main()
