"""Shared utilities: logging, checkpointing, physical constants and conversions."""
import json
import logging
import os
import sys
import time

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_RAW = os.path.join(BASE, "data", "raw")
DATA_PROC = os.path.join(BASE, "data", "processed")
TABLES = os.path.join(BASE, "data", "tables")
LOGS = os.path.join(BASE, "logs")
FIGS = os.path.join(BASE, "output", "figures")
RESULTS = os.path.join(BASE, "output", "results")
for d in (DATA_RAW, DATA_PROC, TABLES, LOGS, FIGS, RESULTS):
    os.makedirs(d, exist_ok=True)

# --- Physical constants (cgs) ---
MSUN_G = 1.98892e33
G_CGS = 6.6743e-8
C_CGS = 2.99792458e10
KPC_CM = 3.0857e21

def r_s_km(m_msun):
    """Schwarzschild radius in km."""
    return 2 * G_CGS * m_msun * MSUN_G / C_CGS**2 / 1e5

def l_edd(m_msun):
    """Eddington luminosity (erg/s) for solar hydrogen composition."""
    return 1.26e38 * m_msun

# --- Crab conversions (approximate monitor calibrations) ---
# MAXI/GSC: Crab ~ 3.7 ph/s/cm2 in 2-20 keV; 2-4 keV ~1.87, 4-10 keV ~1.24, 10-20 keV ~0.29
CRAB_MAXI = {"2-20": 3.7, "2-4": 1.87, "4-10": 1.24, "10-20": 0.29}
# Swift/BAT transient monitor: Crab = 0.220 ct/s/cm2 (15-50 keV)
CRAB_BAT = 0.220
# RXTE/ASM: Crab = 75 ct/s (1.5-12 keV); A(1.5-3)=26.8, B(3-5)=23.3, C(5-12)=25.4
CRAB_ASM = {"sum": 75.0, "a": 26.8, "b": 23.3, "c": 25.4}
# Crab flux in erg/cm2/s: 2-20 keV ~ 2.4e-8; 1.5-12 keV ~ 2.0e-8; 15-50 keV ~ 1.3e-8
CRAB_FLUX = {"maxi_2_20": 2.4e-8, "asm_1p5_12": 2.0e-8, "bat_15_50": 1.3e-8}


def get_logger(name):
    os.makedirs(LOGS, exist_ok=True)
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    fh = logging.FileHandler(os.path.join(LOGS, f"{name}.log"))
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


class Checkpoint:
    """JSON-file based checkpoint: mark atomic units done, skip on rerun."""

    def __init__(self, name):
        self.path = os.path.join(LOGS, f"ckpt_{name}.json")
        self.state = {}
        if os.path.exists(self.path):
            with open(self.path) as f:
                self.state = json.load(f)

    def done(self, key):
        return self.state.get(key, {}).get("status") == "ok"

    def mark(self, key, status="ok", **meta):
        self.state[key] = {"status": status, "t": time.time(), **meta}
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.state, f, indent=1)
        os.replace(tmp, self.path)
