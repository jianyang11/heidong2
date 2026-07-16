"""MAD saturation proxy phi_jet per source from the hard-state Lr-Lx database.

Method (documented in docs/proposal.md sec 3.3):
  P_jet = A * (L_R / 1e30 erg/s)^(12/17)
    flat-spectrum compact-jet scaling L_R ~ P_jet^(17/12) (Blandford & Konigl
    1979; Heinz & Sunyaev 2003; Heinz & Grimm 2005). The normalization A is
    uncertain by ~1 dex; we adopt three calibrations A = {0.4, 2, 8} x 1e36
    erg/s (bracketing Heinz & Grimm 2005 Cyg X-1-anchored and Kording, Fender
    & Migliari 2006 accretion-anchored values) and propagate as systematic.
  Mdot c^2 = L_X,bol / eta ; L_X,bol = 2 * L_X(1-10 keV) (hard-state BC,
    Migliari & Fender 2006). eta = 0.1 for l = L_bol/L_Edd >= 0.01 and
    eta = 0.1 (l/0.01)^0.5 below (radiatively inefficient branch, Narayan &
    Yi 1995 scaling; as in Coriat+11).
  BZ: P_BZ = (kappa/4pi) phi^2 Mdot c^2 x_a^2, x_a = a/(2(1+sqrt(1-a^2))),
    kappa = 0.053 (TNM11). => phi_obs = sqrt(4pi P_jet/(kappa Mdot c^2 x_a^2))
  MAD saturation: phi_max ~= 50.
For sources without spin measurement, a=0.5 assumed (flagged).
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils import BASE, RESULTS, FIGS, TABLES, l_edd, r_s_km, get_logger

log = get_logger("phi_jet")

KAPPA = 0.053
A_CAL = {"low": 0.4e36, "mid": 2e36, "high": 8e36}
BC = 2.0

ALIASES = {
    "GX339-4": "GX 339-4", "GX 339-4": "GX 339-4",
    "MAXI J1820+070": "MAXI J1820+070", "MAXIJ1820+070": "MAXI J1820+070",
    "XTE J1550-564": "XTE J1550-564", "XTEJ1550-564": "XTE J1550-564",
    "GRO J1655-40": "GRO J1655-40", "4U 1543-47": "4U 1543-47",
    "XTE J1859+226": "XTE J1859+226", "GRS 1915+105": "GRS 1915+105",
    "V404 Cyg": "V404 Cyg", "V 404 Cyg": "V404 Cyg",
    "XTE J1118+480": "XTE J1118+480", "XTEJ1118+480": "XTE J1118+480",
    "A0620-00": "A0620-00", "GS 1124-684": "GRS 1124-684",
    "Cyg X-1": "Cyg X-1", "MAXI J1348-630": "MAXI J1348-630",
    "H1743-322": "H1743-322", "H 1743-322": "H1743-322",
    "MAXI J1535-571": "MAXI J1535-571",
    "Swift J1727.8-1613": "Swift J1727.8-1613",
    "MAXI J1659-152": "MAXI J1659-152", "XTE J1752-223": "XTE J1752-223",
}


def x_a(a):
    return a / (2 * (1 + np.sqrt(1 - min(a, 0.998)**2)))


def eta_rad(l_bol_edd):
    return 0.1 if l_bol_edd >= 0.01 else 0.1 * np.sqrt(l_bol_edd / 0.01)


def phi_obs(lr, lx, m, a, A):
    p_jet = A * (lr / 1e30)**(12 / 17)
    l_bol = BC * lx
    l = l_bol / l_edd(m)
    mdot_c2 = l_bol / eta_rad(l)
    return np.sqrt(4 * np.pi * p_jet / (KAPPA * mdot_c2 * x_a(a)**2))


def main():
    sample = pd.read_csv(os.path.join(TABLES, "bh_sample.csv")).set_index("name")
    db = pd.read_csv(os.path.join(
        BASE, "data", "XRB-LrLx_pub", "data", "lrlx_data_BHs.csv"))
    db["name_std"] = db["Name"].map(lambda n: ALIASES.get(str(n).strip()))
    db = db[db["name_std"].notna() & db["uplim"].isna()]
    recs = []
    for name, g in db.groupby("name_std"):
        if name not in sample.index:
            continue
        row = sample.loc[name]
        a = row["spin"] if np.isfinite(row["spin"]) else 0.5
        for _, r in g.iterrows():
            lr, lx = float(r["Lr"]), float(r["Lx"])
            if lr <= 0 or lx <= 0:
                continue
            rec = {"name": name, "Lr": lr, "Lx": lx,
                   "ledd": BC * lx / l_edd(row["M_msun"]),
                   "M": row["M_msun"], "M_err": row["M_err"],
                   "Rs_km": r_s_km(row["M_msun"]),
                   "spin": a, "spin_meas": np.isfinite(row["spin"]),
                   "tier": row["tier"]}
            for k, A in A_CAL.items():
                rec[f"phi_{k}"] = phi_obs(lr, lx, row["M_msun"], a, A)
            recs.append(rec)
    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(RESULTS, "phi_jet_points.csv"), index=False)
    # per-source hard-state maximum phi (use points with ledd in [1e-4, 0.05]
    # to avoid quiescence (different physics) and transition contamination)
    m = (df.ledd > 1e-4) & (df.ledd < 0.05)
    agg = df[m].groupby("name").apply(
        lambda g: pd.Series({
            "phi_max_mid": g["phi_mid"].max(),
            "phi_med_mid": g["phi_mid"].median(),
            "phi_max_low": g["phi_low"].max(),
            "phi_max_high": g["phi_high"].max(),
            "n_pts": len(g),
            "M": g["M"].iloc[0], "M_err": g["M_err"].iloc[0],
            "Rs_km": g["Rs_km"].iloc[0], "spin": g["spin"].iloc[0],
            "spin_meas": bool(g["spin_meas"].iloc[0]),
            "tier": g["tier"].iloc[0],
        }), include_groups=False).reset_index()
    agg.to_csv(os.path.join(RESULTS, "phi_jet_per_source.csv"), index=False)
    log.info("phi per source:\n%s",
             agg[["name", "n_pts", "phi_max_mid", "Rs_km"]].to_string())
    # figures: phi vs ledd trajectories
    fig, ax = plt.subplots(figsize=(8, 6))
    for name, g in df.groupby("name"):
        g = g.sort_values("ledd")
        ax.plot(g["ledd"], g["phi_mid"], "o-", ms=3, lw=0.8, alpha=0.7,
                label=name)
    ax.axhline(50, color="k", ls="--", lw=1.5)
    ax.text(2e-9, 55, r"$\phi_{\rm MAD}\approx 50$ (TNM11)", fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$L_{\rm bol}/L_{\rm Edd}$")
    ax.set_ylabel(r"$\phi_{\rm jet}$ (dimensionless magnetic flux)")
    ax.legend(fontsize=7, ncol=2)
    ax.set_title("Jet-inferred magnetic flux vs Eddington ratio (mid calib.)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "phi_vs_ledd.png"), dpi=140)
    plt.close(fig)
    log.info("done")


if __name__ == "__main__":
    main()
