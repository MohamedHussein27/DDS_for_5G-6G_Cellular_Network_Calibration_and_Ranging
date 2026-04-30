"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_fixed.py

    Description:
        This module provides the top-level wrapper for the bit-true IFFT, applying 
        twiddle conjugation to the core FFT engine. 
        (Plotting utilities have been removed for simulation performance).
"""

import numpy as np

# Import from your split files
from fft_fixed import *

def radix2_dif_ifft_fixed(x_ints: np.ndarray, WL: int, FL: int) -> np.ndarray:
    """
    100% bit-true IFFT wrapper.
    Exact Python translation of radix2_dif_ifft_fixed.m.
    """
    return radix2_dif_fft_fixed(x_ints, WL, FL, is_ifft=True)

# ══════════════════════════════════════════════════════════════════════════════
# Demo
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    N  = 4096
    WL = 16
    FL = 5 

    print("=" * 60)
    print("  Bit-True IFFT Golden Model — Demo (Headless)")
    print(f"  N={N}  WL={WL}  FL={FL}  (Q{WL-FL-1}.{FL})")
    print("=" * 60)

    # ── Build a test input: single tone at bin k=10 ───────────────────────
    scale    = 1 << FL           
    X_float  = np.zeros(N, dtype=np.complex128)
    X_float[10]   = 0.5 + 0j     
    X_float[N-10] = 0.5 + 0j    
    X_in = (X_float * scale).astype(np.int64)  

    print(f"\n── Input: single tone at k=10  (Q11.5 integers, peak={int(0.5*scale)}) ──")
    print(f"  Input Re range : [{np.real(X_in).min()}, {np.real(X_in).max()}]")
    print(f"  Input Im range : [{np.imag(X_in).min()}, {np.imag(X_in).max()}]")

    # ── Run the bit-true IFFT ─────────────────────────────────────────────
    print("\n── Running radix2_dif_ifft_fixed ──")
    x_out = radix2_dif_ifft_fixed(X_in.astype(np.complex128), WL, FL)

    print(f"  Input  length : {len(X_in)}")
    print(f"  Output length : {len(x_out)}")
    print(f"  Output Re range : [{np.real(x_out).min():.0f}, {np.real(x_out).max():.0f}]")
    print(f"  Output Im range : [{np.imag(x_out).min():.0f}, {np.imag(x_out).max():.0f}]")

    print("\n── Execution Complete. ──")