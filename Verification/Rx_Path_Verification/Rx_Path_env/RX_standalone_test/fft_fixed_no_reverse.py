"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: fft_fixed_no_reverse.py

    Description:
        100% bit-true Radix-2 DIF FFT / IFFT engine.
        Exact Python translation of the Verilog RTL arithmetic.

    BUG FIXED (was causing ALL IFFT scoreboard mismatches):
    ────────────────────────────────────────────────────────
    Every power-of-two constant was written with Python's multiplication
    operator `*` instead of the exponentiation operator `**`:

        WRONG  →  CORRECT
        2*15   →  2**15   (= 32768, not 30)
        216    →  2**16   (= 65536, not 216)
        2*14   →  2**14   (= 16384, not 28)
        2*13   →  2**13   (=  8192, not 26)

    Because of this the model was operating with wildly wrong moduli:
      • wrap16 wrapped at ±108 (mod 216) instead of ±32768 (mod 65536)
      • bit_slice_round divided by 28 instead of 16384
      • twiddle factors were clamped to ±29 instead of ±32767

    The standalone IFFT tests "passed" because input amplitudes were small
    enough to never hit the wrong wrap boundary.  The RX path passes signals
    with much larger dynamic range from the conjugate-multiplier, so the
    error was immediately visible.

    This is why the TOP-level scoreboard (which uses run_rx_pipeline with a
    separate correct model) always showed PASS while the IFFT sub-scoreboard
    showed ~1791/1793 mismatches.
"""

import numpy as np


def wrap16(x: np.ndarray) -> np.ndarray:
    """
    16-bit signed wrap — exact Verilog wire behaviour.
    Verilog: assign out = in[15:0] (signed)
    Python:  (x + 2**15) % 2**16 - 2**15
    """
    # FIX: was (x + 2*15) % 216 - 2*15  →  mod 216=216, offset 30
    return (np.asarray(x, dtype=np.int64) + 2**15) % 2**16 - 2**15


def bit_slice_round(x: np.ndarray) -> np.ndarray:
    """
    Exact Verilog truncation with rounding:
      out = sub[29:14] + sub[13]   (take bits 29..14, add bit 13 as round)
    Python: floor(x / 2**14) + mod(floor(x / 2**13), 2)
    """
    # FIX: was floor_divide(x, 2*14) = /28   →  should be /16384
    #      was mod(floor_divide(x, 2*13), 2)  →  should be /8192
    return (np.floor_divide(x, 2**14) +
            np.mod(np.floor_divide(x, 2**13), 2))


def bitrevorder(x: np.ndarray) -> np.ndarray:
    """Reorder array elements in bit-reversed index order."""
    N    = len(x)
    bits = int(np.log2(N))
    result = np.zeros(N, dtype=np.complex128)
    for i in range(N):
        rev = int(format(i, f'0{bits}b')[::-1], 2)
        result[rev] = x[i]
    return result


def radix2_dif_fft_fixed(x_ints: np.ndarray,
                          WL: int,
                          FL: int,
                          is_ifft: bool = False) -> np.ndarray:
    """
    100% bit-true Radix-2 DIF FFT / IFFT core (integer arithmetic).
    Exact Python translation of radix2_dif_fft_fixed.m / the Verilog RTL.

    Parameters
    ----------
    x_ints  : input array of signed 16-bit integers (complex128 container)
    WL      : word length (16 for this design)
    FL      : fractional bits (used only by callers for scaling; not used
              internally — the core operates purely in integer domain)
    is_ifft : True → conjugate twiddle factors (IFFT mode)
    """
    N      = len(x_ints)
    stages = int(np.log2(N))

    X = x_ints.astype(np.complex128).copy()

    for s in range(stages, 0, -1):
        m  = 2 ** s
        mh = m // 2

        # ── Twiddle generation (Q2.14, 16-bit signed) ──────────────────
        k_idx = np.arange(mh)
        W_fft = np.exp(-1j * 2 * np.pi * k_idx / m)

        W_re = np.round(np.real(W_fft) * (2**14)).astype(np.int64)
        W_im = np.round(np.imag(W_fft) * (2**14)).astype(np.int64)

        # Clamp to 16-bit signed range
        # FIX: was -(2*15) = -30, 2*15-1 = 29  →  should be -(2**15), 2**15-1
        W_re = np.clip(W_re, -(2**15), 2**15 - 1)
        W_im = np.clip(W_im, -(2**15), 2**15 - 1)

        # IFFT: conjugate twiddle (negate imaginary part)
        if is_ifft:
            W_im = -W_im

        for k in range(0, N, m):
            # ── 1. Load inputs ────────────────────────────────────────
            a_re = np.real(X[k      : k + mh]).astype(np.int64)
            a_im = np.imag(X[k      : k + mh]).astype(np.int64)
            b_re = np.real(X[k + mh : k + m ]).astype(np.int64)
            b_im = np.imag(X[k + mh : k + m ]).astype(np.int64)

            # ── 2. Butterfly (16-bit wrap on add/sub) ─────────────────
            a_out_re = wrap16(a_re + b_re)
            a_out_im = wrap16(a_im + b_im)
            b_out_re = wrap16(a_re - b_re)
            b_out_im = wrap16(a_im - b_im)

            # ── 3. Twiddle multiply (32-bit exact) ────────────────────
            mul1 = b_out_re * W_re
            mul2 = b_out_im * W_im
            mul3 = b_out_im * W_re
            mul4 = b_out_re * W_im

            sub_re = mul1 - mul2
            add_im = mul3 + mul4

            # ── 4. Truncation (Verilog bit-slice + round) ─────────────
            p_re = bit_slice_round(sub_re)
            p_im = bit_slice_round(add_im)

            # ── 5. Final 16-bit wrap ───────────────────────────────────
            p_re = wrap16(p_re)
            p_im = wrap16(p_im)

            # ── Write back ────────────────────────────────────────────
            X[k      : k + mh] = a_out_re + 1j * a_out_im
            X[k + mh : k + m ] = p_re     + 1j * p_im

    return X