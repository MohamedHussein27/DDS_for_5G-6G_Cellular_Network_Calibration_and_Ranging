import numpy as np
from mytypes import quantize_fixed, reinterpretcast
from fft_fixed import *


# ---------------------------------------------------------------------------
# radix22_dif_ifft_fixed
# ---------------------------------------------------------------------------

def radix22_dif_ifft_fixed(x: np.ndarray, WL: int, FL: int):
    """
    Fixed-point IFFT via circular folding property.

    Mirrors radix22_dif_ifft_fixed.m:
        IFFT(X) ≡ FFT( fold(X) ) / N

    The fold operation is: keep index 0, reverse indices 1..N-1.

    Parameters
    ----------
    x   : Input array (complex). Length must be a power of 2.
    WL  : Word length
    FL  : Fraction length

    Returns
    -------
    X            : IFFT output (reinterpreted bits, see MATLAB logic)
    WL_output    : Output word length
    FL_output    : Output fraction length
    """
    N = len(x)
    integer_part = WL - FL
    growth_bits  = int(np.log2(N))

    # ------------------------------------------------------------------
    # Circular folding: [x[0], x[N-1], x[N-2], ..., x[1]]
    # ------------------------------------------------------------------
    X_folded = x.copy()
    if N > 1:
        X_folded[1:] = x[N - 1:0:-1]   # reverse indices 1..N-1

    # ------------------------------------------------------------------
    # FFT of the folded sequence
    # ------------------------------------------------------------------
    X_fi = radix22_dif_fft_fixed(X_folded, WL, FL)

    # ------------------------------------------------------------------
    # Reinterpret cast (matches MATLAB reinterpretcast logic)
    # ------------------------------------------------------------------
    if integer_part > growth_bits:
        # out_type = numerictype(1, WL, FL + growth_bits)
        dst_FL = FL + growth_bits
        X = reinterpretcast(X_fi,
                            src_WL=WL, src_FL=FL,
                            dst_WL=WL, dst_FL=dst_FL)
        WL_output = WL
        FL_output = WL - 1

    else:
        # out_type2 = numerictype(1, WL, FL + growth_bits)
        dst_FL = FL + growth_bits
        X = reinterpretcast(X_fi,
                            src_WL=WL, src_FL=FL,
                            dst_WL=WL, dst_FL=dst_FL)
        WL_output = WL
        FL_output = FL + growth_bits    # not explicitly set in MATLAB else branch

    return X, WL_output, FL_output