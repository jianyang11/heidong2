"""Download MAXI/GSC and Swift/BAT monitoring light curves for the BH sample.

Checkpointed; concurrency <= 4; logs progress in real time.
Usage: python3 download_xray.py [--pilot]
"""
import os
import sys
import concurrent.futures as cf

import pandas as pd
import requests

from utils import DATA_RAW, TABLES, get_logger, Checkpoint

log = get_logger("download_xray")
ck = Checkpoint("download_xray")

MAXI_URL = "http://maxi.riken.jp/star_data/{sid}/{sid}_g_lc_1day_all.dat"
BAT_URL = "https://swift.gsfc.nasa.gov/results/transients/weak/{name}.lc.txt"
BAT_URL2 = "https://swift.gsfc.nasa.gov/results/transients/{name}.lc.txt"

PILOT = ["GX 339-4", "MAXI J1820+070", "XTE J1550-564"]

ASM_URL = ("https://heasarc.gsfc.nasa.gov/FTP/xte/data/archive/ASMProducts/"
           "definitive_1dwell/lightcurves/xa_{aid}_d1.lc")
ASM_COL_URL = ("https://heasarc.gsfc.nasa.gov/FTP/xte/data/archive/ASMProducts/"
               "definitive_1dwell/colors/xa_{aid}_d1.col")
ASM_IDS = {
    "GX 339-4": "gx339-4", "XTE J1550-564": "xtej1550-564",
    "GRO J1655-40": "groj1655-40", "XTE J1859+226": "xtej1859+226",
    "XTE J1118+480": "xtej1118+480", "4U 1543-47": "x1543-475",
    "GRS 1124-684": "gs1124-684", "V404 Cyg": "gs2023+338",
    "H1743-322": "h1743-322", "GRS 1915+105": "grs1915+105",
    "Cyg X-1": "cygx1", "A0620-00": "x0620-003",
}


def fetch(url, dest, binary=False):
    r = requests.get(url, timeout=300)
    if r.status_code != 200 or len(r.content) < 200:
        return False, f"http {r.status_code} len {len(r.content)}"
    with open(dest, "wb") as f:
        f.write(r.content)
    return True, f"{len(r.content)} bytes"


def do_source(row):
    name = row["name"]
    safe = name.replace(" ", "_")
    results = []
    if isinstance(row.get("maxi_id"), str) and row["maxi_id"]:
        key = f"maxi:{safe}"
        if not ck.done(key):
            dest = os.path.join(DATA_RAW, f"maxi_{safe}.dat")
            ok, msg = fetch(MAXI_URL.format(sid=row["maxi_id"]), dest)
            ck.mark(key, "ok" if ok else "fail", msg=msg)
            log.info("MAXI %s: %s %s", name, "OK" if ok else "FAIL", msg)
            results.append(("maxi", ok))
    if isinstance(row.get("bat_name"), str) and row["bat_name"]:
        key = f"bat:{safe}"
        if not ck.done(key):
            dest = os.path.join(DATA_RAW, f"bat_{safe}.dat")
            ok, msg = fetch(BAT_URL2.format(name=row["bat_name"]), dest)
            if not ok:
                ok, msg = fetch(BAT_URL.format(name=row["bat_name"]), dest)
            ck.mark(key, "ok" if ok else "fail", msg=msg)
            log.info("BAT %s: %s %s", name, "OK" if ok else "FAIL", msg)
            results.append(("bat", ok))
    if name in ASM_IDS:
        key = f"asm:{safe}"
        if not ck.done(key):
            dest = os.path.join(DATA_RAW, f"asm_{safe}.lc")
            ok, msg = fetch(ASM_URL.format(aid=ASM_IDS[name]), dest, binary=True)
            ck.mark(key, "ok" if ok else "fail", msg=msg)
            log.info("ASM %s: %s %s", name, "OK" if ok else "FAIL", msg)
            results.append(("asm", ok))
        key = f"asmcol:{safe}"
        if not ck.done(key):
            dest = os.path.join(DATA_RAW, f"asmcol_{safe}.col")
            ok, msg = fetch(ASM_COL_URL.format(aid=ASM_IDS[name]), dest, binary=True)
            ck.mark(key, "ok" if ok else "fail", msg=msg)
            log.info("ASMCOL %s: %s %s", name, "OK" if ok else "FAIL", msg)
    return name, results


def main():
    df = pd.read_csv(os.path.join(TABLES, "bh_sample.csv"))
    if "--pilot" in sys.argv:
        df = df[df["name"].isin(PILOT)]
    log.info("downloading for %d sources", len(df))
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        for name, res in ex.map(do_source, [r for _, r in df.iterrows()]):
            pass
    log.info("done")


if __name__ == "__main__":
    main()
