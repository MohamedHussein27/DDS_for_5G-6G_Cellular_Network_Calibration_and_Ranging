"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: mixed_ofdm_radar_float_seq.py

    Description:
        TC-RX-005-FLOAT: Mixed OFDM+Radar — Floating Point Verification.
        Bypasses fixed-point quantization to demonstrate zero imaginary leakage.
"""

import numpy as np
import cocotb
from pyuvm import *
from top_seq_item import *
from cocotb.triggers import *

# ─────────────────────────────────────────────────────────────────────────────
# Constants (Kept identical for sizing, but FL is ignored for data generation)
# ─────────────────────────────────────────────────────────────────────────────
N_FFT         = 4096
N_HALF        = 2048
N_OFDM_SC     = 1500
RESET_CYCLES  = 8
DRAIN_CYCLES  = N_FFT + N_HALF * 3 + 256
TARGET_MAX_RX = 225.0 / 256.0  # Floating-point equivalent of the integer scale
AMP_REF_FLOAT = 120.0 / 256.0  # Floating-point equivalent of the reference amplitude

def _build_mixed_rx_frame_float(seed: int = 0) -> list:
    """
    Build a 4096-sample time-domain frame using pure floats (NO QUANTIZATION).
    """
    rng  = np.random.default_rng(seed)
    freq = np.zeros(N_FFT, dtype=np.complex128)

    # OFDM in bins 0..2047
    sc_start   = (N_HALF - N_OFDM_SC) // 2
    sc_end     = sc_start + N_OFDM_SC
    qam_levels = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float64)
    
    freq[sc_start:sc_end] = (1.0 / 7.0) * (
        rng.choice(qam_levels, size=N_OFDM_SC) +
        1j * rng.choice(qam_levels, size=N_OFDM_SC)
    )

    # Chirp in bins 2048..4095
    for k in range(N_HALF, N_FFT):
        k_local = k - N_HALF
        phase   = np.pi * k_local * k_local / N_HALF
        freq[k] = 1.0 * np.exp(1j * phase)

    # IFFT to Time Domain
    td = np.fft.ifft(freq) * N_FFT
    
    # Scale utilizing pure float coefficients
    max_td_val = np.max(np.abs(td))
    scale = TARGET_MAX_RX / max_td_val
    td_scaled = td * scale

    # Return pure floats directly instead of passing through _f2q
    return [(float(s.real), float(s.imag)) for s in td_scaled]


def _build_chirp_reference_float(seed: int = 0) -> list:
    """
    Build 2048 reference entries using pure floats (NO QUANTIZATION).
    """
    ref = []
    for k in range(N_HALF):
        phase = np.pi * k * k / N_HALF
        c     = AMP_REF_FLOAT * np.exp(1j * phase)
        ref.append((float(c.real), float(c.imag)))
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


class mixed_ofdm_radar_float_seq(uvm_sequence):
    """
    Sequence driving unquantized floating-point values to isolate tracking behavior.
    Note: Standard integer RTL will clip/truncate if directly driven with floats;
    this serves as the ideal DSP golden model reference.
    """
    def __init__(self, name="mixed_ofdm_radar_float_seq"):
        super().__init__(name)
        self.num_frames = 1
        self.seed       = 0

    async def body(self):
        for frame_idx in range(self.num_frames):
            cocotb.log.info(f"TC-RX-005-FLOAT frame {frame_idx + 1}: Starting Golden Float Sequence.")

            rx_samples  = _build_mixed_rx_frame_float(seed=self.seed + frame_idx)
            ref_samples = _build_chirp_reference_float(seed=self.seed + frame_idx)

            # Phase 0: Reset
            for cyc in range(RESET_CYCLES):
                await _send(self, f"rst_{frame_idx}_{cyc}", rst_n=0)
            await _send(self, f"rst_release_{frame_idx}", rst_n=1)

            # Phase 1: Reference RAM write
            for addr in range(N_HALF):
                re, im = ref_samples[addr]
                await _send(self, f"ref_wr_{frame_idx}_{addr}", ref_wr=1, ref_wr_re=re, ref_wr_im=im)
            await _send(self, f"ref_wr_done_{frame_idx}")

            # Phase 2: Drive mixed rx_in frame
            for cycle in range(N_FFT):
                re, im = rx_samples[cycle]
                await _send(self, f"rx_{frame_idx}_{cycle}", rx_valid=1, rx_re=re, rx_im=im, ref_wr=0)
            await _send(self, f"rx_done_{frame_idx}")

            # Phase 3: Drain
            for cyc in range(DRAIN_CYCLES):
                await _send(self, f"drain_{frame_idx}_{cyc}")