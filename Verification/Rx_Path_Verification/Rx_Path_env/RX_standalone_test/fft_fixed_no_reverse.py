"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: fft_fixed.py

    Description:
        This module contains the core 100% bit-true Radix-2 DIF FFT engine.
        It precisely mimics the integer arithmetic of the Verilog RTL.
"""

import numpy as np

def wrap16(x: np.ndarray) -> np.ndarray:
    """
    16-bit WRAP (EXACT Verilog wire behavior).
    Translates: mod(x + 2^15, 2^16) - 2^15
    """
    return (np.asarray(x, dtype=np.int64) + 2*15) % 216 - 2*15

def bit_slice_round(x: np.ndarray) -> np.ndarray:
    """
    EXACT VERILOG TRUNCATION MODEL.
    Translates: floor(x / 2^14) + mod(floor(x / 2^13), 2)
    sub_re[29:14] + sub_re[13]
    """
    return np.floor_divide(x, 2*14) + np.mod(np.floor_divide(x, 2*13), 2)

def bitrevorder(x: np.ndarray) -> np.ndarray:
    """
    Reorder array elements in bit-reversed index order.
    """
    N = len(x)
    bits = int(np.log2(N))
    result = np.zeros(N, dtype=np.complex128)
    for i in range(N):
        rev = int(format(i, f'0{bits}b')[::-1], 2)
        result[rev] = x[i]
    return result

def radix2_dif_fft_fixed(x_ints: np.ndarray, WL: int, FL: int, is_ifft: bool = False) -> np.ndarray:
    """
    100% bit-true FFT / IFFT core (integer arithmetic).
    Exact Python translation of radix2_dif_fft_fixed.m.
    """
    N = len(x_ints)
    stages = int(np.log2(N))

    # Work in 64-bit float/complex internally for safe 32-bit wire simulation
    X = x_ints.astype(np.complex128).copy()

    for s in range(stages, 0, -1):
        m  = 2**s
        mh = m // 2

        # ===============================
        # Twiddle generation (Q2.14)
        # ===============================
        k_idx = np.arange(mh)
        W_fft = np.exp(-1j * 2 * np.pi * k_idx / m)

        W_re = np.round(np.real(W_fft) * (2**14)).astype(np.int64)
        W_im = np.round(np.imag(W_fft) * (2**14)).astype(np.int64)

        # Clamp to 16-bit signed
        W_re = np.clip(W_re, -(2*15), 2*15 - 1)
        W_im = np.clip(W_im, -(2*15), 2*15 - 1)

        # Twiddle Selection (FFT vs IFFT)
        if is_ifft:
            W_im = -W_im

        for k in range(0, N, m):
            # =========================================================
            # 1. LOAD INPUTS
            # =========================================================
            a_re = np.real(X[k : k+mh]).astype(np.int64)
            a_im = np.imag(X[k : k+mh]).astype(np.int64)
            b_re = np.real(X[k+mh : k+m]).astype(np.int64)
            b_im = np.imag(X[k+mh : k+m]).astype(np.int64)

            # =========================================================
            # 2. BUTTERFLY (16-bit WRAP ONLY)
            # =========================================================
            a_out_re = wrap16(a_re + b_re)
            a_out_im = wrap16(a_im + b_im)
            b_out_re = wrap16(a_re - b_re)
            b_out_im = wrap16(a_im - b_im)

            # =========================================================
            # 3. MULTIPLIER (32-bit EXACT like Verilog)
            # =========================================================
            mul1 = b_out_re * W_re
            mul2 = b_out_im * W_im
            mul3 = b_out_im * W_re
            mul4 = b_out_re * W_im

            sub_re = mul1 - mul2
            add_im = mul3 + mul4

            # =========================================================
            # 4. TRUNCATION EXACT VERILOG STYLE
            # =========================================================
            p_re = bit_slice_round(sub_re)
            p_im = bit_slice_round(add_im)

            # =========================================================
            # 5. FINAL WRAP TO 16-bit (overflow allowed)
            # =========================================================
            p_re = wrap16(p_re)
            p_im = wrap16(p_im)

            # =========================================================
            # WRITE BACK
            # =========================================================
            X[k : k+mh]   = a_out_re + 1j * a_out_im
            X[k+mh : k+m] = p_re + 1j * p_im

    return X