"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: radar_range_profile_seq.py

    Description:
        TC-RX-007: Radar Range Profile — Loaded Reference.
        Feature: F-RX-01, F-RX-02

    STIMULUS:
        Pre-load the 2048-entry reference RAM with a TX chirp reference
        derived from the golden MATLAB/Python model.
        Apply a matching received chirp signal as rx_in.
        Capture radar_out.

    EXPECTED / PASS:
        Bit-exact match between radar_out and the MATLAB range-profile
        golden reference (Error = 0). All 2048 samples verified.

    Sequence structure
    ──────────────────────────────────────────────────────────────
    Phase 0 — Reset                (RESET_CYCLES)
    Phase 1 — Reference RAM write  (2048 cycles, TX chirp reference)
    Phase 2 — Chirp rx_in frame    (4096 cycles)
    Phase 3 — Drain                (DRAIN_CYCLES)
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
MAX_VAL      =  (1 << (WL - 1)) - 1
MIN_VAL      = -(1 << (WL - 1))
RESET_CYCLES = 8
DRAIN_CYCLES = N_FFT + N_HALF * 3 + 256


def _f2q(f: float) -> int:
    return max(MIN_VAL, min(MAX_VAL, int(round(f * (1 << FL)))))


def _build_chirp_rx_frame() -> list:
    """
    Build a 4096-sample time-domain frame with chirp energy in
    radar bins 2048..4095 only. OFDM bins 0..2047 = 0.
    """
    freq = np.zeros(N_FFT, dtype=np.complex128)
    amp  = MAX_VAL / 8.0
    for k in range(N_HALF, N_FFT):
        k_local = k - N_HALF
        phase   = np.pi * k_local * k_local / N_HALF
        freq[k] = amp * np.exp(1j * phase)
    td = np.fft.ifft(freq) * N_FFT
    return [(_f2q(float(s.real)), _f2q(float(s.imag))) for s in td]


def _build_tx_chirp_reference() -> list:
    """
    Build the 2048-entry TX chirp reference that matches the chirp
    spectrum used in the rx frame (same phase law, same amplitude).
    Returns list of (re_int, im_int) for reference RAM addresses 0..2047.
    """
    amp = MAX_VAL / 8.0
    ref = []
    for k in range(N_HALF):
        phase = np.pi * k * k / N_HALF
        c     = amp * np.exp(1j * phase)
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


# ─────────────────────────────────────────────────────────────────────────────
# TC-RX-007 sequence
# ─────────────────────────────────────────────────────────────────────────────
class radar_range_profile_seq(uvm_sequence):
    """
    TC-RX-007: Radar Range Profile — Loaded Reference.
    Loads a TX chirp reference into the RAM, then drives a matching
    received chirp. Scoreboard checks bit-exact match on all 2048
    radar_out samples against the golden range-profile model.

    Attributes:
        num_frames : int — number of frames (default 1)
    """

    def __init__(self, name="radar_range_profile_seq"):
        super().__init__(name)
        self.num_frames = 1

    async def body(self):

        chirp_samples = _build_chirp_rx_frame()
        ref_samples   = _build_tx_chirp_reference()

        for frame_idx in range(self.num_frames):

            cocotb.log.info(
                f"TC-RX-007 frame {frame_idx + 1}/{self.num_frames}: "
                f"starting — Radar range profile, loaded reference."
            )

            # ── Phase 0: Reset ────────────────────────────────────────────
            for cyc in range(RESET_CYCLES):
                await _send(self, f"rst_{frame_idx}_{cyc}", rst_n=0)
            await _send(self, f"rst_release_{frame_idx}", rst_n=1)

            cocotb.log.info(
                f"TC-RX-007 frame {frame_idx + 1}: reset complete."
            )

            # ── Phase 1: Reference RAM write ──────────────────────────────
            cocotb.log.info(
                f"TC-RX-007 frame {frame_idx + 1}: "
                f"writing {N_HALF} TX chirp reference entries (ref_wr_en=1)..."
            )

            for addr in range(N_HALF):
                re, im = ref_samples[addr]
                await _send(
                    self,
                    f"ref_wr_{frame_idx}_{addr}",
                    ref_wr=1,
                    ref_re=re,
                    ref_im=im,
                )

            await _send(self, f"ref_wr_done_{frame_idx}")

            cocotb.log.info(
                f"TC-RX-007 frame {frame_idx + 1}: reference RAM write complete."
            )

            # ── Phase 2: Drive chirp rx_in frame ──────────────────────────
            cocotb.log.info(
                f"TC-RX-007 frame {frame_idx + 1}: "
                f"driving {N_FFT} chirp samples (rx_valid_in=1)..."
            )

            for cycle in range(N_FFT):
                re, im = chirp_samples[cycle]
                await _send(
                    self,
                    f"chirp_{frame_idx}_{cycle}",
                    rx_valid=1,
                    rx_re=re,
                    rx_im=im,
                    ref_wr=0,
                )

            await _send(self, f"chirp_done_{frame_idx}")

            cocotb.log.info(
                f"TC-RX-007 frame {frame_idx + 1}: "
                f"chirp frame complete. Entering drain phase..."
            )

            # ── Phase 3: Drain ────────────────────────────────────────────
            for cyc in range(DRAIN_CYCLES):
                await _send(self, f"drain_{frame_idx}_{cyc}")

            cocotb.log.info(
                f"TC-RX-007 frame {frame_idx + 1}: drain complete — frame done."
            )
