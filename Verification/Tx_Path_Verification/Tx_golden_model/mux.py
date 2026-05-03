"""
mux.py
======
Python equivalent of mux.m

Combines positive-frequency and negative-frequency streams into a
standard FFT-order spectrum buffer (ready to feed directly into IFFT).

Inputs
------
stream_pos : Positive-frequency data  (0 → +Fs/2), length N/2
stream_neg : Negative-frequency data  (-Fs/2 → 0), length N/2
N          : Total FFT size

Outputs
-------
X_out      : Combined spectrum in standard FFT order, shape (N,)
valid_out  : Boolean "data valid" strobe, shape (N,)
"""

import numpy as np


def mux(stream_pos: np.ndarray,
        stream_neg: np.ndarray,
        N: int):
    """
    Hardware MUX: combine positive and negative frequency streams.

    Mirrors the hardware pipeline in mux.m exactly:
      - Phase 1 (k = 0 .. N/2-1) : write from stream_pos
      - Phase 2 (k = N/2 .. N-1) : write from stream_neg

    Parameters
    ----------
    stream_pos : ndarray, complex, length N/2
    stream_neg : ndarray, complex, length N/2
    N          : int, total FFT size (must be even)

    Returns
    -------
    X_out     : complex ndarray, shape (N,)
    valid_out : bool ndarray, shape (N,)
    """
    X_out     = np.zeros(N, dtype=complex)
    valid_out = np.zeros(N, dtype=bool)

    ptr_pos = 0   # Python uses 0-based indexing  (MATLAB ptr_pos started at 1)
    ptr_neg = 0

    for k in range(N):
        if k < N // 2:
            # Phase 1: Positive Frequencies
            X_out[k] = stream_pos[ptr_pos]
            ptr_pos += 1
        else:
            # Phase 2: Negative Frequencies
            X_out[k] = stream_neg[ptr_neg]
            ptr_neg += 1

        valid_out[k] = True

    return X_out, valid_out
