"""Interpolated cross-correlation (ICCF, Gaskell & Peterson 1987) with
FR/RSS Monte-Carlo errors (Peterson et al. 1998). Positive lag => series-2
lags series-1 (radio lags X-ray if x=X-ray, y=radio).
"""
import numpy as np


def iccf(t1, f1, t2, f2, lags):
    """Two-way interpolation CCF averaged (standard ICCF)."""
    r = np.full(len(lags), np.nan)
    for k, lag in enumerate(lags):
        cc = []
        # interpolate series 2 at t1 + lag
        m = (t1 + lag >= t2.min()) & (t1 + lag <= t2.max())
        if m.sum() > 5:
            y2 = np.interp(t1[m] + lag, t2, f2)
            a, b = f1[m], y2
            if a.std() > 0 and b.std() > 0:
                cc.append(np.corrcoef(a, b)[0, 1])
        # interpolate series 1 at t2 - lag
        m = (t2 - lag >= t1.min()) & (t2 - lag <= t1.max())
        if m.sum() > 5:
            y1 = np.interp(t2[m] - lag, t1, f1)
            a, b = y1, f2[m]
            if a.std() > 0 and b.std() > 0:
                cc.append(np.corrcoef(a, b)[0, 1])
        if cc:
            r[k] = np.mean(cc)
    return r


def centroid_lag(lags, r, frac=0.8):
    if np.all(~np.isfinite(r)):
        return np.nan, np.nan
    imax = np.nanargmax(r)
    rmax = r[imax]
    m = r >= frac * rmax
    # keep contiguous region around peak
    idx = np.where(m)[0]
    grp = idx[np.searchsorted(idx, imax) == np.arange(len(idx))] if False else None
    lo = imax
    while lo > 0 and m[lo - 1]:
        lo -= 1
    hi = imax
    while hi < len(m) - 1 and m[hi + 1]:
        hi += 1
    sel = slice(lo, hi + 1)
    cen = np.nansum(lags[sel] * r[sel]) / np.nansum(r[sel])
    return cen, rmax


def frrss(t1, f1, e1, t2, f2, e2, lags, n=500, seed=0):
    """FR/RSS MC: returns centroid distribution."""
    rng = np.random.default_rng(seed)
    cens = []
    for _ in range(n):
        i1 = np.sort(rng.integers(0, len(t1), len(t1)))
        i1 = np.unique(i1)
        i2 = np.sort(rng.integers(0, len(t2), len(t2)))
        i2 = np.unique(i2)
        ff1 = f1[i1] + rng.normal(0, e1[i1])
        ff2 = f2[i2] + rng.normal(0, e2[i2])
        r = iccf(t1[i1], ff1, t2[i2], ff2, lags)
        c, rm = centroid_lag(lags, r)
        if np.isfinite(c) and rm > 0.3:
            cens.append(c)
    return np.array(cens)


def measure_lag(t1, f1, e1, t2, f2, e2, lag_min=-30, lag_max=30, dlag=0.5,
                nmc=400, seed=0):
    lags = np.arange(lag_min, lag_max + dlag, dlag)
    r = iccf(np.asarray(t1, float), np.asarray(f1, float),
             np.asarray(t2, float), np.asarray(f2, float), lags)
    cen, rmax = centroid_lag(lags, r)
    cens = frrss(np.asarray(t1, float), np.asarray(f1, float),
                 np.asarray(e1, float), np.asarray(t2, float),
                 np.asarray(f2, float), np.asarray(e2, float), lags,
                 n=nmc, seed=seed)
    if len(cens) < nmc * 0.3:
        return {"lag": cen, "rmax": rmax, "err_lo": np.nan, "err_hi": np.nan,
                "lags": lags, "r": r, "cens": cens, "quality": "C"}
    lo, med, hi = np.percentile(cens, [15.87, 50, 84.13])
    qual = "A" if (len(cens) > nmc * 0.7 and rmax > 0.6) else "B"
    return {"lag": med, "rmax": rmax, "err_lo": med - lo, "err_hi": hi - med,
            "lags": lags, "r": r, "cens": cens, "quality": qual}


if __name__ == "__main__":
    # unit test: recover a known lag
    rng = np.random.default_rng(1)
    t = np.arange(0, 200, 1.0)
    sig = np.exp(-0.5 * ((t - 80) / 25)**2)
    t2 = np.arange(0, 200, 3.0)
    true_lag = 7.0
    sig2 = np.interp(t2 - true_lag, t, sig)
    f1 = sig + rng.normal(0, 0.03, len(t))
    f2 = sig2 + rng.normal(0, 0.05, len(t2))
    res = measure_lag(t, f1, np.full(len(t), .03), t2, f2, np.full(len(t2), .05))
    print(f"true lag 7.0 -> measured {res['lag']:.2f} "
          f"+{res['err_hi']:.2f}/-{res['err_lo']:.2f} rmax={res['rmax']:.2f} "
          f"quality={res['quality']}")
    assert abs(res["lag"] - true_lag) < 1.5, "FAILED to recover lag"
    print("unit test PASSED")
