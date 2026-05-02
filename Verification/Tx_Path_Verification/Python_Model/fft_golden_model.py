"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: ifft_golden_model.py
 
    Description:
        Exact Python translation of radix2_dif_fft_fixed.m (parameterized version).
 
        Mirrors every line of the MATLAB code:
            - TW_FL = 14  (Q2.14 twiddles, hardcoded same as MATLAB)
            - SCALE_TW   = 2^TW_FL   = 16384
            - HALF_SHIFT = 2^(TW_FL-1) = 8192
            - wrap_nbit() wraps to any WL signed integer (same as Verilog wire)
            - Butterfly: parameterized wrap on add/sub outputs
            - Multiplier: exact 32-bit products, nearest-integer truncation
            - Output: bit-reversal permutation
"""
 
import numpy as np
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Helper — parameterized N-bit signed wrap
# Exact translation of MATLAB local function wrap_nbit(in, WL)
# ─────────────────────────────────────────────────────────────────────────────
def wrap_nbit(x: np.ndarray, WL: int) -> np.ndarray:
    """
    Wrap integer array to signed WL-bit range [-2^(WL-1), 2^(WL-1)-1].
 
    MATLAB equivalent:
        half_range = 2^(WL-1);
        full_range = 2^WL;
        out = mod(in + half_range, full_range) - half_range;
    """
    half_range = 1 << (WL - 1)
    full_range = 1 << WL
    x = np.asarray(x, dtype=np.int64)
    return (x + half_range) % full_range - half_range
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Helper — bit-reversal permutation
# Equivalent to MATLAB's bitrevorder()
# ─────────────────────────────────────────────────────────────────────────────
def _bitrevorder(x: np.ndarray) -> np.ndarray:
    N    = len(x)
    bits = int(np.log2(N))
    out  = np.zeros(N, dtype=np.complex128)
    for i in range(N):
        rev      = int(f'{i:0{bits}b}'[::-1], 2)
        out[rev] = x[i]
    return out
 
 
# ─────────────────────────────────────────────────────────────────────────────
# radix2_dif_fft_fixed
# Exact translation of radix2_dif_fft_fixed.m (parameterized version)
# ─────────────────────────────────────────────────────────────────────────────
def radix2_dif_fft_fixed(x_ints: np.ndarray,
                          WL: int,
                          is_ifft: bool = False) -> np.ndarray:
    """
    Parameterized bit-true Radix-2 DIF FFT / IFFT core.
    Exact Python translation of radix2_dif_fft_fixed.m.
 
    Parameters
    ──────────
    x_ints  : complex integer input array, length N (must be power of 2)
    WL      : data word length (e.g. 16)
    is_ifft : False → forward FFT (default)
              True  → IFFT  (conjugates twiddle imaginary parts)
 
    Returns
    ───────
    X_ints  : complex integer output array in bit-reversed order
    """
 
    # TW_FL hardcoded to 14 — same as MATLAB "TW_FL = 14; % Default to Q2.14"
    TW_FL      = 14
    SCALE_TW   = 1 << TW_FL         # 2^TW_FL   = 16384
    HALF_SHIFT = 1 << (TW_FL - 1)   # 2^(TW_FL-1) = 8192
 
    N      = len(x_ints)
    stages = int(np.log2(N))
    assert 2 ** stages == N, "Input length must be a power of 2."
 
    # Work with float64 internally — same as MATLAB "X = double(x_ints)"
    X = x_ints.astype(np.complex128).copy()
 
    for s in range(stages, 0, -1):      # s = stages:-1:1
        m  = 1 << s                     # 2^s
        mh = m >> 1                     # m/2
 
        # ── 1. Generate base FFT twiddles ────────────────────────────────
        # W_fft = exp(-1j*2*pi*(0:mh-1)/m).'
        k_idx    = np.arange(mh)
        W_fft    = np.exp(-1j * 2 * np.pi * k_idx / m)
 
        # Scale twiddles dynamically based on TW_FL
        W_re_raw = np.round(np.real(W_fft) * SCALE_TW).astype(np.int64)
        W_im_raw = np.round(np.imag(W_fft) * SCALE_TW).astype(np.int64)
 
        # ── 2. Twiddle selection — FFT vs IFFT ───────────────────────────
        # W_re = W_re_raw  (unchanged)
        # if is_ifft: W_im = -W_im_raw  else: W_im = W_im_raw
        W_re = W_re_raw
        W_im = -W_im_raw if is_ifft else W_im_raw
 
        # ── 3. Butterfly groups ──────────────────────────────────────────
        # for k = 1:m:N  (MATLAB 1-indexed → Python 0-indexed: range(0, N, m))
        for k in range(0, N, m):
            # Butterfly inputs
            a_re = np.real(X[k      : k + mh]).astype(np.int64)
            a_im = np.imag(X[k      : k + mh]).astype(np.int64)
            b_re = np.real(X[k + mh : k + m ]).astype(np.int64)
            b_im = np.imag(X[k + mh : k + m ]).astype(np.int64)
 
            # ── Butterfly stage (parameterized wrapping) ─────────────────
            a_out_re = wrap_nbit(a_re + b_re, WL)
            a_out_im = wrap_nbit(a_im + b_im, WL)
            b_out_re = wrap_nbit(a_re - b_re, WL)
            b_out_im = wrap_nbit(a_im - b_im, WL)
 
            # ── Multiplier stage (exact math) ────────────────────────────
            mul1   = b_out_re * W_re
            mul2   = b_out_im * W_im
            sub_re = mul1 - mul2        # re_out = (b_re * W_re) - (b_im * W_im)
 
            mul3   = b_out_im * W_re
            mul4   = b_out_re * W_im
            add_im = mul3 + mul4        # im_out = (b_im * W_re) + (b_re * W_im)
 
            # ── Fixed-point truncation ────────────────────────────────────
            # p_re = floor(sub_re / SCALE_TW) + mod(floor(sub_re / HALF_SHIFT), 2)
            # p_im = floor(add_im / SCALE_TW) + mod(floor(add_im / HALF_SHIFT), 2)
            p_re = sub_re // SCALE_TW + (sub_re // HALF_SHIFT) % 2
            p_im = add_im // SCALE_TW + (add_im // HALF_SHIFT) % 2
 
            # Wrap multiplier outputs back to parameterized WL format
            p_re = wrap_nbit(p_re, WL)
            p_im = wrap_nbit(p_im, WL)
 
            # ── Write back ───────────────────────────────────────────────
            X[k      : k + mh] = a_out_re + 1j * a_out_im
            X[k + mh : k + m ] = p_re     + 1j * p_im
 
    # Bit-reversed output — same as MATLAB bitrevorder()
    #return _bitrevorder(X)
    return X
 
 
# ─────────────────────────────────────────────────────────────────────────────
# radix2_dif_ifft_fixed — bit-true IFFT wrapper
# Exact translation of radix2_dif_ifft_fixed.m
# ─────────────────────────────────────────────────────────────────────────────
def radix2_dif_ifft_fixed(x_ints: np.ndarray, WL: int) -> np.ndarray:
    """
    Bit-true IFFT wrapper.
    Calls radix2_dif_fft_fixed with is_ifft=True — same as the MATLAB wrapper
    that passes is_ifft=1 to apply twiddle conjugation exactly as the RTL does.
 
    Parameters
    ──────────
    x_ints : complex integer frequency-domain input, length N (power of 2)
    WL     : data word length (e.g. 16)
 
    Returns
    ───────
    X_ints : complex integer time-domain output in bit-reversed order
    """
    return radix2_dif_fft_fixed(x_ints, WL, is_ifft=True)