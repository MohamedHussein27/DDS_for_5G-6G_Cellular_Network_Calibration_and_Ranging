"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: mixed_ofdm_radar_small_seq.py

    Description:
        TC-RX-005-SMALL: Micro-Signal Verification.
        Drives the exact mixed OFDM+Radar frame, but the maximum amplitudes 
        are massively constrained to absolutely guarantee no internal FFT overflow.
"""

import numpy as np
import cocotb
from pyuvm import *
from top_seq_item import *
from cocotb.triggers import *

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
WL           = 16
FL           = 8
N_FFT        = 4096
N_HALF       = 2048
N_OFDM_SC    = 1500
MAX_VAL      =  (1 << (WL - 1)) - 1
MIN_VAL      = -(1 << (WL - 1))
RESET_CYCLES = 8
DRAIN_CYCLES = N_FFT + N_HALF * 3 + 256

# =============================================================================
# ULTRA-SMALL SCALING CONSTANTS (NO OVERFLOW GUARANTEE)
# =============================================================================
# Reduced from 225 down to 10.
# The maximum integer pushing into your Q8.8 16-bit wires will be +/- 10.
TARGET_MAX_RX_INT = 10

# Reduced from 120.0 down to 10.0.
AMP_REF_FLOAT = 10.0


def _f2q(f: float) -> int:
    return max(MIN_VAL, min(MAX_VAL, int(round(f * (1 << FL)))))


def _build_mixed_rx_frame_small(seed: int = 0) -> list:
    rng  = np.random.default_rng(seed)
    freq = np.zeros(N_FFT, dtype=np.complex128)

    sc_start   = (N_HALF - N_OFDM_SC) // 2
    sc_end     = sc_start + N_OFDM_SC
    qam_levels = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float64)
    
    freq[sc_start:sc_end] = (1.0 / 7.0) * (
        rng.choice(qam_levels, size=N_OFDM_SC) +
        1j * rng.choice(qam_levels, size=N_OFDM_SC)
    )

    for k in range(N_HALF, N_FFT):
        k_local = k - N_HALF
        phase   = np.pi * k_local * k_local / N_HALF
        freq[k] = 1.0 * np.exp(1j * phase)

    td = np.fft.ifft(freq) * N_FFT
    
    max_td_val = np.max(np.abs(td))
    scale = (TARGET_MAX_RX_INT / float(1 << FL)) / max_td_val
    td_scaled = td * scale

    return [(_f2q(float(s.real)), _f2q(float(s.imag))) for s in td_scaled]


def _build_chirp_reference_small(seed: int = 0) -> list:
    ref = []
    for k in range(N_HALF):
        phase = np.pi * k * k / N_HALF
        c     = AMP_REF_FLOAT * np.exp(1j * phase)
        ref.append((_f2q(float(c.real)), _f2q(float(c.imag))))
    return ref


async def _send(seq, name,
                rst_n=1, rx_valid=0, rx_re=0, rx_im=0,
                ref_wr=0, ref_re=0, ref_im=0):
    req = rx_item.create(name)
    await seq.start_item(req)
    req.rst_n       = rst_n
    req.rx_valid_in = rx_valid
    req.rx_in_re    = rx_re
    req.rx_in_im    = rx_im
    req.ref_wr_en   = ref_wr
    req.ref_wr_re   = ref_re
    req.ref_wr_im   = ref_im
    await seq.finish_item(req)


class mixed_ofdm_radar_small_seq(uvm_sequence):
    def __init__(self, name="mixed_ofdm_radar_small_seq"):
        super().__init__(name)
        self.num_frames = 1
        self.seed       = 0

    async def body(self):
        for frame_idx in range(self.num_frames):
            cocotb.log.info(
                f"TC-RX-005-SMALL frame {frame_idx + 1}: "
                f"starting — MICRO-SIGNAL Mode, seed={self.seed + frame_idx}."
            )

            rx_samples  = _build_mixed_rx_frame_small(seed=self.seed + frame_idx)
            ref_samples = _build_chirp_reference_small(seed=self.seed + frame_idx)

            for cyc in range(RESET_CYCLES):
                await _send(self, f"rst_{frame_idx}_{cyc}", rst_n=0)
            await _send(self, f"rst_release_{frame_idx}", rst_n=1)

            for addr in range(N_HALF):
                re, im = ref_samples[addr]
                await _send(self, f"ref_wr_{frame_idx}_{addr}", ref_wr=1, ref_re=re, ref_im=im)
            await _send(self, f"ref_wr_done_{frame_idx}")

            for cycle in range(N_FFT):
                re, im = rx_samples[cycle]
                await _send(self, f"rx_{frame_idx}_{cycle}", rx_valid=1, rx_re=re, rx_im=im)
            await _send(self, f"rx_done_{frame_idx}")

            for cyc in range(DRAIN_CYCLES):
                await _send(self, f"drain_{frame_idx}_{cyc}")