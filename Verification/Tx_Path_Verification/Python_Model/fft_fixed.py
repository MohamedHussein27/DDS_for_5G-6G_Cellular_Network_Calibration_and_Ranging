"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: fft_fixed.py

    Description:
        This module contains the core 100% bit-true Radix-2^2 DIF FFT engine.
        It precisely mimics the integer arithmetic of the Verilog RTL, including 
        Q2.14 twiddle factor generation, 16-bit signed wrapping on butterfly 
        add/sub stages, 32-bit multiplier truncation with nearest-integer rounding,
        and final bit-reversal permutation.
"""

import numpy as np

def _wrap_16bit(x: np.ndarray) -> np.ndarray:
    """
    Wrap any integer array to the range [-32768, 32767] exactly like a
    16-bit Verilog wire (modular arithmetic, no saturation).
    """
    return (np.asarray(x, dtype=np.int64) + 32768) % 65536 - 32768

def _bitrevorder(x: np.ndarray) -> np.ndarray:
    """
    Reorder array elements in bit-reversed index order.
    """
    N      = len(x)
    bits   = int(np.log2(N))
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
    N      = len(x_ints)
    stages = int(np.log2(N))
    assert 2 ** stages == N, "Input length must be a power of 2."

    # Work in 64-bit float to safely simulate 32-bit wire arithmetic
    X = x_ints.astype(np.complex128).copy()

    for s in range(stages, 0, -1):
        m  = 2 ** s
        mh = m // 2

        # Generate base FFT twiddles in Q2.14
        k_idx    = np.arange(mh)
        W_fft    = np.exp(-1j * 2 * np.pi * k_idx / m)

        W_re_raw = np.round(np.real(W_fft) * 16384).astype(np.int64)
        W_im_raw = np.round(np.imag(W_fft) * 16384).astype(np.int64)

        W_re_raw = np.clip(W_re_raw, -32768, 32767)
        W_im_raw = np.clip(W_im_raw, -32768, 32767)

        # Twiddle selection — conjugate for IFFT
        W_re = W_re_raw
        W_im = -W_im_raw if is_ifft else W_im_raw

        for k in range(0, N, m):
            a_re = np.real(X[k      : k + mh]).astype(np.int64)
            a_im = np.imag(X[k      : k + mh]).astype(np.int64)
            b_re = np.real(X[k + mh : k + m ]).astype(np.int64)
            b_im = np.imag(X[k + mh : k + m ]).astype(np.int64)

            # Butterfly stage (pure add / sub, 16-bit wrap)
            a_out_re = _wrap_16bit(a_re + b_re)
            a_out_im = _wrap_16bit(a_im + b_im)
            b_out_re = _wrap_16bit(a_re - b_re)
            b_out_im = _wrap_16bit(a_im - b_im)

            # Multiplier stage (32-bit exact)
            mul1   = b_out_re * W_re
            mul2   = b_out_im * W_im
            sub_re = mul1 - mul2

            mul3   = b_out_im * W_re
            mul4   = b_out_re * W_im
            add_im = mul3 + mul4

            # Verilog truncation: sub_re[29:14] + sub_re[13]
            p_re = sub_re // 16384 + (sub_re // 8192) % 2
            p_im = add_im // 16384 + (add_im // 8192) % 2

            # Wrap multiplier outputs back to 16-bit signed
            p_re = _wrap_16bit(p_re)
            p_im = _wrap_16bit(p_im)

            X[k      : k + mh] = a_out_re + 1j * a_out_im
            X[k + mh : k + m ] = p_re     + 1j * p_im

    X_ints = _bitrevorder(X)
    return X_ints