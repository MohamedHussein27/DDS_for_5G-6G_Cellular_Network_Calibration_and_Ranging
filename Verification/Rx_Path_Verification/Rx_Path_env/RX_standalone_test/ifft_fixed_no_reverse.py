"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_fixed_no_reverse.py

    Description:
        This module provides the top-level wrapper for the bit-true IFFT.
"""

import numpy as np
from fft_fixed_no_reverse import radix2_dif_fft_fixed

def radix2_dif_ifft_fixed(x_ints: np.ndarray, WL: int, FL: int, is_ifft: bool = True) -> np.ndarray:
    """
    100% BIT-TRUE IFFT WRAPPER
    Instead of circular folding, this passes an 'is_ifft = True' flag to the 
    FFT core, instructing it to use Twiddle Conjugation (exactly like the RTL).
    """
    is_ifft = is_ifft
    X_ints = radix2_dif_fft_fixed(x_ints, WL, FL, is_ifft)
    
    return X_ints