"""
mytypes.py
==========
Python equivalent of mytypes.m

Defines fixed-point quantization helpers for the DSP golden model.
Since Python has no native fi() type, we emulate MATLAB's fi() behavior
using integer arithmetic (numpy int32/int64) with explicit word-length (WL)
and fraction-length (FL) parameters.

Usage:
    T = mytypes('fixed', dds_frac_bits=16)
    T = mytypes('double')
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Fixed-point helper functions (replaces fi / fimath in MATLAB)
# ---------------------------------------------------------------------------

def quantize_fixed(x: np.ndarray, WL: int, FL: int,
                   signed: bool = True,
                   rounding: str = 'floor',
                   overflow: str = 'wrap') -> np.ndarray:
    """
    Quantize a floating-point array to a fixed-point representation.

    Parameters
    ----------
    x        : Input array (float / complex)
    WL       : Total word length (bits)
    FL       : Fractional bits
    signed   : True for signed (2's complement), False for unsigned
    rounding : 'floor'     -> truncation  (matches MATLAB Floor)
               'convergent'-> round-half-to-even (matches MATLAB Convergent)
               'nearest'   -> round-half-away-from-zero (matches MATLAB Nearest)
    overflow : 'wrap'      -> modular wrap  (matches MATLAB Wrap)
               'saturate'  -> clamp to min/max (matches MATLAB Saturate)

    Returns
    -------
    Quantized values as a float array (storing the fixed-point *value*,
    not the raw integer).  For raw integers multiply by 2**FL.
    """
    is_complex = np.iscomplexobj(x)

    def _quantize_real(v):
        scale = 2.0 ** FL

        # --- Rounding ---
        if rounding == 'floor':
            v_int = np.floor(v * scale)
        elif rounding == 'convergent':
            # Round-half-to-even (banker's rounding)
            v_int = np.round(v * scale)          # numpy uses banker's rounding
        elif rounding == 'nearest':
            v_int = np.floor(v * scale + 0.5)
        else:
            raise ValueError(f"Unknown rounding method: {rounding}")

        # --- Overflow ---
        if signed:
            max_int =  2 ** (WL - 1) - 1
            min_int = -2 ** (WL - 1)
        else:
            max_int = 2 ** WL - 1
            min_int = 0

        if overflow == 'saturate':
            v_int = np.clip(v_int, min_int, max_int)
        elif overflow == 'wrap':
            total = 2 ** WL
            v_int = v_int % total
            if signed:
                v_int[v_int > max_int] -= total
        else:
            raise ValueError(f"Unknown overflow mode: {overflow}")

        return v_int / scale

    if is_complex:
        return _quantize_real(np.real(x)) + 1j * _quantize_real(np.imag(x))
    else:
        return _quantize_real(x)


def reinterpretcast(x: np.ndarray,
                    src_WL: int, src_FL: int,
                    dst_WL: int, dst_FL: int) -> np.ndarray:
    """
    Reinterpret the raw integer bits of x under a new (WL, FL) type,
    i.e. the integer bit pattern is preserved and only the scaling changes.

    Equivalent to MATLAB's reinterpretcast(x, out_type).
    """
    if np.iscomplexobj(x):
        real_part = np.round(np.real(x) * (2.0 ** src_FL)).astype(np.int64)
        imag_part = np.round(np.imag(x) * (2.0 ** src_FL)).astype(np.int64)
        scale = 2.0 ** dst_FL
        return real_part / scale + 1j * (imag_part / scale)
    else:
        raw = np.round(x * (2.0 ** src_FL)).astype(np.int64)
        return raw / (2.0 ** dst_FL)


# ---------------------------------------------------------------------------
# Type descriptor dataclass
# ---------------------------------------------------------------------------

@dataclass
class TypeDescriptor:
    """Holds WL/FL parameters for each signal in the DDS/FFT chain."""
    dtype: str          # 'double', 'single', or 'fixed'

    # DDS path
    M_WL:        int = 32
    M_FL:        int = 0
    acc_WL:      int = 32
    acc_FL:      int = 0
    phase_WL:    int = 32
    phase_FL:    int = 0
    lut_addr_WL: int = 16
    lut_addr_FL: int = 0

    # Sine LUT / DDS output  (Q1.fracBits)
    sine_WL:     int = 18   # fracBits + 2 (default fracBits=16)
    sine_FL:     int = 16
    dds_out_WL:  int = 18
    dds_out_FL:  int = 16

    # fimath equivalents
    rounding: str = 'floor'
    overflow: str = 'wrap'


def mytypes(dt: str, dds_frac_bits: Optional[int] = None) -> TypeDescriptor:
    """
    Python equivalent of mytypes.m

    Parameters
    ----------
    dt            : 'double' | 'single' | 'fixed'
    dds_frac_bits : Fractional bits for the sine LUT (fixed mode only).
                    Defaults to 16 if None.
    """
    if dt in ('double', 'single'):
        # In double/single mode all types are floating-point;
        # WL/FL are irrelevant for quantization but kept for interface compat.
        return TypeDescriptor(dtype=dt)

    elif dt == 'fixed':
        frac_bits = dds_frac_bits if dds_frac_bits is not None else 16
        WL = frac_bits + 2          # e.g. 18 for frac_bits=16

        return TypeDescriptor(
            dtype='fixed',
            # Phase path (unsigned 32-bit, no fraction)
            M_WL=32,        M_FL=0,
            acc_WL=32,      acc_FL=0,
            phase_WL=32,    phase_FL=0,
            # LUT address (unsigned 16-bit)
            lut_addr_WL=16, lut_addr_FL=0,
            # Sine amplitude Q1.frac_bits -> WL = frac_bits+2
            sine_WL=WL,     sine_FL=frac_bits,
            dds_out_WL=WL,  dds_out_FL=frac_bits,
            # fimath: Floor rounding, Wrap overflow (matches MATLAB model)
            rounding='floor',
            overflow='wrap',
        )

    else:
        raise ValueError(f"Unknown data type mode: '{dt}'")
