import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# RTL / Q11.5 constants  (shared across all sequences)
# ─────────────────────────────────────────────────────────────────────────────
N       = 4096
WL      = 16
Q_FRAC  = 5
MAX_VAL = (1 << (WL - 1)) - 1   #  32767
MIN_VAL = -(1 << (WL - 1))      # -32768

def _f2q(f):
    """float → signed Q11.5 integer, saturated to [MIN_VAL, MAX_VAL]."""
    return max(MIN_VAL, min(MAX_VAL, int(round(f * (1 << Q_FRAC)))))

def _frame(arr):
    """numpy complex array → list of (re_q11.5, im_q11.5) signed int pairs."""
    return [(_f2q(c.real), _f2q(c.imag)) for c in arr]