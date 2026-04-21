"""
dds_core.py
===========
Python equivalent of dds_core.m

Hardware-accurate Direct Digital Synthesizer (DDS) model.

Inputs
------
M        : Tuning word array  (integer frequency control words)
Nacc     : Phase accumulator bit-width
LUT_bits : Number of bits used to address the sine LUT (must be <= Nacc)
DT_Mode  : 'double' | 'single' | 'fixed'
frac_bits: (optional) Fractional bits for sine LUT in fixed mode (default 16)

Output
------
dds_out  : Sine-wave samples (float array, length = len(M))
"""

import numpy as np
from mytypes import mytypes, quantize_fixed


def dds_core(M: np.ndarray,
             Nacc: int,
             LUT_bits: int,
             DT_Mode: str,
             frac_bits: int = 16) -> np.ndarray:
    """
    Hardware-accurate DDS golden model.

    Mirrors dds_core.m exactly, including:
      - uint64 phase accumulator with natural modular wrap
      - Phase truncation to LUT_bits MSBs
      - Fixed-point sine LUT quantization (when DT_Mode='fixed')
    """
    # ------------------------------------------------------------------
    # 1. Load types
    # ------------------------------------------------------------------
    T = mytypes(DT_Mode, dds_frac_bits=frac_bits)

    assert LUT_bits <= Nacc, "LUT_bits must be <= Nacc"

    # ------------------------------------------------------------------
    # 2. Cast inputs
    # ------------------------------------------------------------------
    M_cast = np.asarray(M, dtype=np.uint64)
    Ns     = len(M_cast)

    # ------------------------------------------------------------------
    # 3. Sine LUT  (Phase-to-Amplitude Converter)
    # ------------------------------------------------------------------
    LUT_size = int(2 ** LUT_bits)
    addr_ROM = np.arange(LUT_size, dtype=np.float64)
    sine_LUT_float = np.sin(2.0 * np.pi * addr_ROM / LUT_size)

    if DT_Mode == 'fixed':
        # Quantize to Q(1, frac_bits) – matches MATLAB fi(…, T.sine_LUT)
        sine_LUT = quantize_fixed(
            sine_LUT_float,
            WL=T.sine_WL,
            FL=T.sine_FL,
            signed=True,
            rounding=T.rounding,   # 'floor'
            overflow=T.overflow,   # 'wrap'
        )
    else:
        # 'double' or 'single' – no quantization
        sine_LUT = sine_LUT_float.astype(
            np.float32 if DT_Mode == 'single' else np.float64
        )

    # ------------------------------------------------------------------
    # 4. Phase Accumulator
    #    Uses uint64 arithmetic so wrap is automatic (matches MATLAB uint64)
    # ------------------------------------------------------------------
    acc_mask = np.uint64(2 ** Nacc - 1)   # Nacc-bit mask

    acc   = np.uint64(0)
    phase = np.zeros(Ns, dtype=np.uint64)

    for n in range(Ns):
        phase[n] = acc
        acc = np.uint64(acc + M_cast[n]) & acc_mask

    # ------------------------------------------------------------------
    # 5. Phase Truncation  (keep only the top LUT_bits)
    # ------------------------------------------------------------------
    addr_shift = Nacc - LUT_bits
    lut_addr   = (phase >> np.uint64(addr_shift)).astype(np.int64)

    # ------------------------------------------------------------------
    # 6. ROM Lookup
    # ------------------------------------------------------------------
    dds_out = sine_LUT[lut_addr]

    return dds_out
