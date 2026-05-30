"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: null_input_seq.py

    Description:
        TC-RX-001: Null Input — All Zeros.
        Feature: F-RX-05

    STIMULUS:
        Set rx_in_re = 0, rx_in_im = 0 for all 4096 samples.
        Assert rx_valid_in for one full frame.
        No reference RAM write (ref_wr_en stays 0).

    EXPECTED / PASS:
        All 2048 ofdm_out samples = 0.
        All 2048 radar_out samples = 0.
        Both valid counts correct (exactly 2048 pulses each).

    Sequence structure
    ──────────────────────────────────────────────────────────────
    Phase 0 — Reset        (RESET_CYCLES)
    Phase 1 — Null frame   (4096 cycles, rx_valid_in=1, re=0, im=0)
    Phase 2 — Drain        (DRAIN_CYCLES idle cycles)
"""

import cocotb
from pyuvm import *
from top_seq_item import *
from cocotb.triggers import *

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
N_FFT        = 4096
N_HALF       = 2048
RESET_CYCLES = 8
DRAIN_CYCLES = N_FFT + N_HALF * 3 + 256


# ─────────────────────────────────────────────────────────────────────────────
# Helper — send one rx_item
# ─────────────────────────────────────────────────────────────────────────────
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
# TC-RX-001 sequence
# ─────────────────────────────────────────────────────────────────────────────
class null_input_seq(uvm_sequence):
    """
    TC-RX-001: Null Input — All Zeros.
    Drives rx_in = 0+0j for one full 4096-sample frame with no reference
    RAM write. Expects all ofdm_out and radar_out samples to be zero.

    Attributes:
        num_frames : int — number of frames to run (default 1)
    """

    def __init__(self, name="null_input_seq"):
        super().__init__(name)
        self.num_frames = 1

    async def body(self):

        for frame_idx in range(self.num_frames):

            cocotb.log.info(
                f"TC-RX-001 frame {frame_idx + 1}/{self.num_frames}: "
                f"starting — Null input, all zeros."
            )

            # ── Phase 0: Reset ────────────────────────────────────────────
            for cyc in range(RESET_CYCLES):
                await _send(self, f"rst_{frame_idx}_{cyc}", rst_n=0)
            await _send(self, f"rst_release_{frame_idx}", rst_n=1)

            cocotb.log.info(
                f"TC-RX-001 frame {frame_idx + 1}: reset complete."
            )

            # ── Phase 1: Null rx_in frame (4096 zero samples) ─────────────
            cocotb.log.info(
                f"TC-RX-001 frame {frame_idx + 1}: "
                f"driving {N_FFT} zero samples (rx_valid_in=1)..."
            )

            for cycle in range(N_FFT):
                await _send(
                    self,
                    f"null_{frame_idx}_{cycle}",
                    rx_valid=1,
                    rx_re=0,
                    rx_im=0,
                    ref_wr=0,
                )

            await _send(self, f"null_done_{frame_idx}")

            cocotb.log.info(
                f"TC-RX-001 frame {frame_idx + 1}: "
                f"null frame complete. Entering drain phase..."
            )

            # ── Phase 2: Drain ────────────────────────────────────────────
            for cyc in range(DRAIN_CYCLES):
                await _send(self, f"drain_{frame_idx}_{cyc}")

            cocotb.log.info(
                f"TC-RX-001 frame {frame_idx + 1}: drain complete — frame done."
            )
