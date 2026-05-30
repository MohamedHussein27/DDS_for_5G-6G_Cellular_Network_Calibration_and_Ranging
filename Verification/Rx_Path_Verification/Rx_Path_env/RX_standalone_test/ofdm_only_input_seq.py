"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ofdm_only_input_seq.py

    Description:
        TC-RX-003: OFDM-Only Input (Radar Bins = 0).
        Feature: F-RX-01

    STIMULUS:
        Apply an OFDM signal occupying bins 0–2047 only.
        Bins 2048–4095 are set to zero (no chirp energy).
        Reference RAM is never written (stays zero).

    EXPECTED / PASS:
        ofdm_out matches the golden model (FFT → bit_rev → bins 0..2047).
        radar_out = 0 on all 2048 samples (zero radar bins → zero multiply
        input → zero IFFT → zero output).
        No cross-path contamination.

    Sequence structure
    ──────────────────────────────────────────────────────────────
    Phase 0 — Reset        (RESET_CYCLES)
    Phase 1 — OFDM rx_in   (4096 cycles)
        Frequency frame: OFDM subcarriers in bins 0–2047, zeros in 2048–4095.
        IFFT → time domain → Q8.8.
    Phase 2 — Drain        (DRAIN_CYCLES)
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
N_OFDM_SC    = 1500          # active subcarriers (centred in bins 0..2047)
MAX_VAL      =  (1 << (WL - 1)) - 1
MIN_VAL      = -(1 << (WL - 1))
RESET_CYCLES = 8
DRAIN_CYCLES = N_FFT + N_HALF * 3 + 256


def _f2q(f: float) -> int:
    return max(MIN_VAL, min(MAX_VAL, int(round(f * (1 << FL)))))


def _build_ofdm_rx_frame(seed: int = 0) -> list:
    """
    Build a 4096-sample time-domain frame with OFDM energy only in
    frequency bins 0–2047 (radar bins 2048–4095 stay zero).

    Method:
        1. Randomly assign 256-QAM symbols to N_OFDM_SC active subcarriers
           centred within bins 0..2047. All other bins = 0.
        2. IFFT → time domain → Q8.8.
    """
    rng       = np.random.default_rng(seed)
    freq      = np.zeros(N_FFT, dtype=np.complex128)

    # Centre the active subcarriers within the OFDM half (bins 0..2047)
    sc_start  = (N_HALF - N_OFDM_SC) // 2
    sc_end    = sc_start + N_OFDM_SC

    # 256-QAM: real and imag each drawn from {±1, ±3, ±5, ±7} (normalised)
    qam_levels = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float64)
    re_syms    = rng.choice(qam_levels, size=N_OFDM_SC)
    im_syms    = rng.choice(qam_levels, size=N_OFDM_SC)
    amp        = (MAX_VAL / 8.0) / 7.0   # scale to avoid saturation
    freq[sc_start:sc_end] = amp * (re_syms + 1j * im_syms)

    # Bins 2048..4095 remain 0  ← radar region intentionally empty
    td = np.fft.ifft(freq) * N_FFT

    return [(_f2q(float(s.real)), _f2q(float(s.imag))) for s in td]


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


# ─────────────────────────────────────────────────────────────────────────────
# TC-RX-003 sequence
# ─────────────────────────────────────────────────────────────────────────────
class ofdm_only_input_seq(uvm_sequence):
    """
    TC-RX-003: OFDM-Only Input (Radar Bins = 0).
    Drives an OFDM signal in bins 0–2047 only. Radar bins are zero.
    Reference RAM is never written.

    Attributes:
        num_frames : int — number of frames (default 1)
        seed       : int — RNG seed for OFDM symbol generation (default 0)
    """

    def __init__(self, name="ofdm_only_input_seq"):
        super().__init__(name)
        self.num_frames = 1
        self.seed       = 0

    async def body(self):

        for frame_idx in range(self.num_frames):

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}/{self.num_frames}: "
                f"starting — OFDM-only input (radar bins = 0), seed={self.seed + frame_idx}."
            )

            ofdm_samples = _build_ofdm_rx_frame(seed=self.seed + frame_idx)

            # ── Phase 0: Reset ────────────────────────────────────────────
            for cyc in range(RESET_CYCLES):
                await _send(self, f"rst_{frame_idx}_{cyc}", rst_n=0)
            await _send(self, f"rst_release_{frame_idx}", rst_n=1)

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: reset complete."
            )

            # ── Phase 1: Drive OFDM frame ─────────────────────────────────
            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: "
                f"driving {N_FFT} OFDM samples (rx_valid_in=1)..."
            )

            for cycle in range(N_FFT):
                re, im = ofdm_samples[cycle]
                await _send(
                    self,
                    f"ofdm_{frame_idx}_{cycle}",
                    rx_valid=1,
                    rx_re=re,
                    rx_im=im,
                    ref_wr=0,
                )

            await _send(self, f"ofdm_done_{frame_idx}")

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: "
                f"OFDM frame complete. Entering drain phase..."
            )

            # ── Phase 2: Drain ────────────────────────────────────────────
            for cyc in range(DRAIN_CYCLES):
                await _send(self, f"drain_{frame_idx}_{cyc}")

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: drain complete — frame done."
            )
