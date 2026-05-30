"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: rtl_hex_chirp_seq.py

    Description:
        TC-RX-003: RTL Hex File Stimulus, Chirp Signal In.
        Feature: F-RX-03

    STIMULUS:
        Reference RAM is pre-loaded from:
            rtl_rx_radar_re.hex  (2048 entries, 16-bit signed Q8.8)
            rtl_rx_radar_im.hex  (2048 entries, 16-bit signed Q8.8)
        One full 4096-sample rx_in frame is driven from:
            rtl_rx_in_re.hex     (4096 entries, 16-bit signed Q8.8)
            rtl_rx_in_im.hex     (4096 entries, 16-bit signed Q8.8)
        All hex values are two's-complement 16-bit signed integers.

    EXPECTED / PASS:
        Depends on the RTL golden outputs paired with the hex files.
        The scoreboard / checker compares radar_out and ofdm_out against
        the expected values produced by the RTL simulation that generated
        these stimulus files.

    Sequence structure  (mirrors zero_ref_chirp_seq.py style)
    ──────────────────────────────────────────────────────────────
    Phase 0 — Reset
        Assert rst_n=0 for RESET_CYCLES clock periods, then release.

    Phase 1 — Reference RAM write  (2048 cycles)
        Write all 2048 reference locations from the radar hex files.
        ref_wr_en=1 each cycle, rx_valid_in=0 throughout.

    Phase 2 — rx_in frame  (4096 cycles)
        Drive rx_valid_in=1 with samples from the rx_in hex files.
        ref_wr_en=0 throughout.

    Phase 3 — Pipeline drain  (DRAIN_CYCLES idle cycles)
        All inputs idle.  Allows FFT + bit_rev + IFFT + bit_rev to flush.
"""

import os
import numpy as np
import cocotb
import logging
from pyuvm import *
from top_seq_item import *
import pyuvm
import cocotb
from cocotb.triggers import *

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
WL           = 16
FL           = 8          # Q8.8 fixed-point
N_FFT        = 4096
N_HALF       = 2048
RESET_CYCLES = 8
# Conservative drain: FFT(4096) + bit_rev + IFFT(2048) + bit_rev + margins
DRAIN_CYCLES = N_FFT + N_HALF * 3 + 256

# ─────────────────────────────────────────────────────────────────────────────
# Default hex file paths  (override via HEX_DIR env-var or seq attributes)
# ─────────────────────────────────────────────────────────────────────────────
_DEFAULT_HEX_DIR = os.path.dirname(os.path.abspath(__file__))

RX_IN_RE_HEX   = os.path.join(_DEFAULT_HEX_DIR, "rtl_rx_in_re.hex")
RX_IN_IM_HEX   = os.path.join(_DEFAULT_HEX_DIR, "rtl_rx_in_im.hex")
RADAR_REF_RE_HEX = os.path.join(_DEFAULT_HEX_DIR, "rtl_mux_out_re.hex")
RADAR_REF_IM_HEX = os.path.join(_DEFAULT_HEX_DIR, "rtl_mux_out_im.hex")


# ─────────────────────────────────────────────────────────────────────────────
# Helper — load a hex file into a list of signed 16-bit integers
# ─────────────────────────────────────────────────────────────────────────────
def _load_hex(path: str) -> list:
    """
    Read a file of newline-separated 4-digit hex strings and return a list
    of signed 16-bit (two's complement) integers.

    Accepts both Windows (CRLF) and Unix (LF) line endings.
    Blank lines are silently skipped.
    """
    values = []
    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            unsigned = int(line, 16) & 0xFFFF        # mask to 16 bits
            # Two's complement conversion
            signed   = unsigned if unsigned < 0x8000 else unsigned - 0x10000
            values.append(signed)
    return values


# ─────────────────────────────────────────────────────────────────────────────
# Helper — send one rx_item
# ─────────────────────────────────────────────────────────────────────────────
async def _send(seq, name,
                rst_n    = 1,
                rx_valid = 0, rx_re = 0, rx_im = 0,
                ref_wr   = 0, ref_re = 0, ref_im = 0):
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
class rtl_hex_chirp_seq(uvm_sequence):
    """
    TC-RX-003: RTL Hex File Stimulus, Chirp Signal In.

    Loads all stimulus from the four RTL-generated hex files:
        rtl_rx_in_re.hex     — 4096 rx_in real  samples
        rtl_rx_in_im.hex     — 4096 rx_in imag  samples
        rtl_rx_radar_re.hex  — 2048 reference RAM real  entries
        rtl_rx_radar_im.hex  — 2048 reference RAM imag  entries

    Attributes (configurable before start):
        num_frames      : int — how many times to replay the stimulus (default 1)
        rx_in_re_path   : str — override path for rx_in real  hex file
        rx_in_im_path   : str — override path for rx_in imag  hex file
        ref_re_path     : str — override path for reference real  hex file
        ref_im_path     : str — override path for reference imag  hex file
    """

    def __init__(self, name="rtl_hex_chirp_seq"):
        super().__init__(name)
        self.num_frames    = 1
        self.rx_in_re_path = RX_IN_RE_HEX
        self.rx_in_im_path = RX_IN_IM_HEX
        self.ref_re_path   = RADAR_REF_RE_HEX
        self.ref_im_path   = RADAR_REF_IM_HEX

    async def body(self):

        # ══════════════════════════════════════════════════════════
        # Pre-load all hex files once — shared across all frames
        # ══════════════════════════════════════════════════════════
        cocotb.log.info(
            "TC-RX-003: Loading hex stimulus files..."
        )

        rx_re_samples  = _load_hex(self.rx_in_re_path)
        rx_im_samples  = _load_hex(self.rx_in_im_path)
        ref_re_samples = _load_hex(self.ref_re_path)
        ref_im_samples = _load_hex(self.ref_im_path)

        # Validate lengths
        assert len(rx_re_samples) == N_FFT, (
            f"rx_in_re hex: expected {N_FFT} entries, got {len(rx_re_samples)}"
        )
        assert len(rx_im_samples) == N_FFT, (
            f"rx_in_im hex: expected {N_FFT} entries, got {len(rx_im_samples)}"
        )
        assert len(ref_re_samples) == N_HALF, (
            f"radar_ref_re hex: expected {N_HALF} entries, got {len(ref_re_samples)}"
        )
        assert len(ref_im_samples) == N_HALF, (
            f"radar_ref_im hex: expected {N_HALF} entries, got {len(ref_im_samples)}"
        )

        cocotb.log.info(
            f"TC-RX-003: Hex files loaded — "
            f"{N_FFT} rx_in samples, {N_HALF} reference entries."
        )

        # ══════════════════════════════════════════════════════════
        # Repeat for the configured number of frames
        # ══════════════════════════════════════════════════════════
        for frame_idx in range(self.num_frames):

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}/{self.num_frames}: "
                f"starting — RTL hex stimulus."
            )

            # ══════════════════════════════════════════════════════════
            # PHASE 0: Reset
            # Assert rst_n=0 for RESET_CYCLES, then release.
            # ══════════════════════════════════════════════════════════
            for cyc in range(RESET_CYCLES):
                await _send(self, f"rst_{frame_idx}_{cyc}", rst_n=0)

            # One idle cycle with reset released
            await _send(self, f"rst_release_{frame_idx}", rst_n=1)

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: reset complete."
            )

            # ══════════════════════════════════════════════════════════
            # PHASE 1: Reference RAM write  (2048 cycles)
            # Write all N_HALF reference entries from the radar hex files.
            # ref_wr_en=1, rx_valid_in=0 throughout.
            # ══════════════════════════════════════════════════════════
            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: "
                f"writing {N_HALF} reference RAM entries (ref_wr_en=1)..."
            )

            for addr in range(N_HALF):
                await _send(
                    self,
                    f"ref_wr_{frame_idx}_{addr}",
                    rst_n    = 1,
                    rx_valid = 0,
                    ref_wr   = 1,
                    ref_re   = ref_re_samples[addr],
                    ref_im   = ref_im_samples[addr],
                )

            # One idle cycle to de-assert ref_wr_en
            await _send(self, f"ref_wr_done_{frame_idx}", rst_n=1)

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: reference RAM write complete."
            )

            # ══════════════════════════════════════════════════════════
            # PHASE 2: Drive rx_in frame  (4096 cycles)
            # rx_valid_in=1 with samples from the rx_in hex files.
            # ref_wr_en=0 throughout.
            # ══════════════════════════════════════════════════════════
            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: "
                f"driving {N_FFT} rx_in samples (rx_valid_in=1)..."
            )

            for cycle in range(N_FFT):
                await _send(
                    self,
                    f"rx_{frame_idx}_{cycle}",
                    rst_n    = 1,
                    rx_valid = 1,
                    rx_re    = rx_re_samples[cycle],
                    rx_im    = rx_im_samples[cycle],
                    ref_wr   = 0,
                )

            # De-assert rx_valid_in at end of frame
            await _send(self, f"rx_done_{frame_idx}", rst_n=1)

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: "
                f"rx_in frame complete. Entering drain phase..."
            )

            # ══════════════════════════════════════════════════════════
            # PHASE 3: Pipeline drain
            # All inputs idle.  Wait for both ofdm_valid_out and
            # radar_valid_out to fully arrive from the pipeline.
            # ══════════════════════════════════════════════════════════
            for cyc in range(DRAIN_CYCLES):
                await _send(self, f"drain_{frame_idx}_{cyc}", rst_n=1)

            cocotb.log.info(
                f"TC-RX-003 frame {frame_idx + 1}: "
                f"drain complete — frame done."
            )
