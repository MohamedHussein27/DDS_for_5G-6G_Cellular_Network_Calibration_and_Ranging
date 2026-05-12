"""
mux.py
======
Python equivalent of the isac_mux Verilog RTL.

Combines positive-frequency (OFDM) and negative-frequency (Radar/Chirp) 
streams into a standard FFT-order spectrum buffer.

Mirrors the Verilog State Machine exactly:
  - State 1: 2048 bins from OFDM
  - State 2:  381 bins of Zeros (padding)
  - State 3: 1667 bins from Radar RAM (Shifted right by 7 bits)
"""

import numpy as np

def mux(stream_pos: np.ndarray, stream_neg: np.ndarray, N: int = 4096):
    """
    Hardware MUX: combine positive and negative frequency streams.

    Parameters
    ----------
    stream_pos : ndarray, complex (Expects at least 2048 OFDM samples)
    stream_neg : ndarray, complex (Expects at least 1667 Radar samples)
    N          : int, total FFT size (must be 4096 to match RTL state machine)

    Returns
    -------
    X_out     : complex ndarray, shape (N,)
    valid_out : bool ndarray, shape (N,)
    """
    # Initialize the output array (Defaults to exactly 0, taking care of State 2 naturally)
    X_out     = np.zeros(N, dtype=complex)
    valid_out = np.ones(N, dtype=bool)  # Entire generated frame is valid

    # Ensure inputs are large enough to avoid indexing errors
    if len(stream_pos) < 2048:
        stream_pos = np.pad(stream_pos, (0, 2048 - len(stream_pos)))
    if len(stream_neg) < 1667:
        stream_neg = np.pad(stream_neg, (0, 1667 - len(stream_neg)))

    # =========================================================
    # APPLY HARDWARE SHIFT: FFT Input >> 7
    # =========================================================
    # Extract real and imag, shift as integers, and recombine
    stream_neg_re_shifted = np.real(stream_neg).astype(np.int64) >> 7
    stream_neg_im_shifted = np.imag(stream_neg).astype(np.int64) >> 7
    stream_neg_shifted = stream_neg_re_shifted + 1j * stream_neg_im_shifted

    # =========================================================
    # STATE 1: Output 2048 OFDM bins (rd_cnt == 0 to 2047)
    # =========================================================
    ofdm_len = 2048
    X_out[0 : ofdm_len] = stream_pos[0 : ofdm_len]

    # =========================================================
    # STATE 2: Output 381 Zeros (rd_cnt == 0 to 380)
    # =========================================================
    zero_len = 381
    # X_out is already initialized to zero, so this segment is naturally correct.

    # =========================================================
    # STATE 3: Output 1667 Radar bins (rd_cnt == 0 to 1666)
    # =========================================================
    radar_len = 1667
    radar_start_idx = ofdm_len + zero_len  # 2048 + 381 = 2429
    
    # Reads the first 1667 items of the *shifted* stream
    X_out[radar_start_idx : radar_start_idx + radar_len] = stream_neg_shifted[0 : radar_len]

    return X_out, valid_out