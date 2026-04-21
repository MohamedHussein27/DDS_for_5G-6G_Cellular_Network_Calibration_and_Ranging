"""
fft_fixed.py
============
Python equivalents of:
  - radix22_dif_fft_fixed.m
  - radix22_dif_ifft_fixed.m

Implements a fixed-point Radix-2² DIF FFT with 2 internal guard bits,
matching the active (non-commented) MATLAB implementation.

Key design choices that mirror the MATLAB code:
  - 18-bit internal guard-RAM  : WL_internal = WL + 2, FL_internal = FL + 2
  - Twiddle factors            : Q(1, WL-2), i.e. numerictype(1, WL, WL-2)
  - Rounding                   : Convergent (banker's rounding)
  - Overflow                   : Saturate
  - Final output               : Quantized back to strict (WL, FL)
  - Bit-reversal               : Applied after all butterfly stages
"""

import numpy as np
from mytypes import quantize_fixed, reinterpretcast


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _bit_reverse_order(x: np.ndarray) -> np.ndarray:
    """Return x permuted into bit-reversed order (equivalent to bitrevorder)."""
    N = len(x)
    assert N & (N - 1) == 0, "N must be a power of 2"
    bits = int(np.log2(N))
    indices = np.arange(N, dtype=np.int64)
    rev = np.zeros(N, dtype=np.int64)
    for i in range(bits):
        rev = (rev << 1) | (indices & 1)
        indices >>= 1
    return x[rev]


def _quantize_internal(x: np.ndarray, WL: int, FL: int) -> np.ndarray:
    """Quantize complex array to (WL, FL) with convergent rounding + saturate."""
    return quantize_fixed(x, WL=WL, FL=FL,
                          signed=True,
                          rounding='convergent',
                          overflow='saturate')


# ---------------------------------------------------------------------------
# radix22_dif_fft_fixed
# ---------------------------------------------------------------------------

def radix22_dif_fft_fixed(x: np.ndarray, WL: int, FL: int) -> np.ndarray:
    """
    Fixed-point Radix-2² DIF FFT with 2 guard bits.

    Mirrors the active MATLAB implementation (the "Added Guard bits" version).

    Parameters
    ----------
    x   : Input array (complex or real). Length must be a power of 2.
    WL  : Output word length  (e.g. 16)
    FL  : Output fraction length (e.g. 8)

    Returns
    -------
    X   : FFT output, quantized to (WL, FL), in bit-reversed order.
    """
    N = len(x)
    stages = int(np.log2(N))
    assert 2 ** stages == N, "Input length must be a power of 2"

    # ------------------------------------------------------------------
    # 1. 18-bit internal guard-RAM
    # ------------------------------------------------------------------
    WL_internal = WL + 2
    FL_internal = FL + 2

    # Initialize X in 18-bit guard RAM
    X = _quantize_internal(x.astype(complex), WL_internal, FL_internal)

    # ------------------------------------------------------------------
    # 2. Twiddle factors: Q(1, WL-2) i.e. numerictype(1, WL, WL-2)
    # ------------------------------------------------------------------
    WL_twiddle = WL
    FL_twiddle = WL - 2     # Same as MATLAB T_twiddle = numerictype(1, WL, WL-2)

    # ------------------------------------------------------------------
    # 3. DIF butterfly stages (from stage = stages down to 1)
    # ------------------------------------------------------------------
    for s in range(stages, 0, -1):
        m  = int(2 ** s)
        mh = m // 2

        # Twiddle vector W[k] = exp(-j*2*pi*k/m) for k=0..mh-1
        k_vec  = np.arange(mh, dtype=np.float64)
        W_raw  = np.exp(-1j * 2.0 * np.pi * k_vec / m)

        # Quantize twiddles to (WL, WL-2)
        W = _quantize_internal(W_raw, WL_twiddle, FL_twiddle)

        for k in range(0, N, m):
            a = X[k      : k + mh].copy()
            b = X[k + mh : k + m ].copy()

            # Butterfly: sum and difference, then multiply diff by twiddle
            sum_ab  = a + b
            diff_ab = (a - b) * W

            # Write back to 18-bit guard RAM
            X[k      : k + mh] = _quantize_internal(sum_ab,  WL_internal, FL_internal)
            X[k + mh : k + m ] = _quantize_internal(diff_ab, WL_internal, FL_internal)

    # ------------------------------------------------------------------
    # 4. Bit-reversed output
    # ------------------------------------------------------------------
    X = _bit_reverse_order(X)

    # ------------------------------------------------------------------
    # 5. Squeeze back to strict (WL, FL) output
    # ------------------------------------------------------------------
    X = _quantize_internal(X, WL, FL)

    return X



