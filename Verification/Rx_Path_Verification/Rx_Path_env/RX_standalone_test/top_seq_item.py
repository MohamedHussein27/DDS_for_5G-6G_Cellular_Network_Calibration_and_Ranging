"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: rx_seq_item.py

    Description:
        Transaction-level packet for the RX_TOP block.
        Uses multiple inheritance: UVM item + CRV randomization.

        Word format: WL=16 signed fixed-point for all data fields.

    Review notes vs. original top_seq_item.py
    ──────────────────────────────────────────
    1. Renamed class rx_item → rx_item  (no change, kept consistent)
    2. ref_wr_re / ref_wr_im: removed from add_rand (correct — these are
       driven manually by sequences, not randomized, to control the exact
       reference loaded into the RAM).
    3. Added convert2string_output() for scoreboard-side debug logging
       (mirrors convert2string_stimulus which exists for driver-side).
    4. All other fields and constraints are correct as-is.
"""

import cocotb
from cocotb_coverage.crv import *
from pyuvm import *
import random

class rx_item(uvm_sequence_item, Randomized):

    def __init__(self, name="rx_item"):
        uvm_sequence_item.__init__(self, name)
        Randomized.__init__(self)

        # =========================================================
        # INPUT VARIABLES  (driven by the driver)
        # =========================================================
        self.rst_n       = 1

        # Channel input path
        self.rx_valid_in = 0
        self.rx_in_re    = 0
        self.rx_in_im    = 0

        # Reference RAM write path
        self.ref_wr_en   = 0
        self.ref_wr_re   = 0    # Not randomized — set manually by sequences
        self.ref_wr_im   = 0    # Not randomized — set manually by sequences

        # =========================================================
        # OUTPUT VARIABLES  (captured by the monitor, not randomized)
        # =========================================================
        # Communication (OFDM) output
        self.ofdm_valid_out  = 0
        self.ofdm_out_re     = 0
        self.ofdm_out_im     = 0

        # Radar (range-profile) output
        self.radar_valid_out = 0
        self.radar_out_re    = 0
        self.radar_out_im    = 0

        # =========================================================
        # CRV SETUP
        # =========================================================
        self.add_rand("rst_n",       [0, 1])
        self.add_rand("rx_valid_in", [0, 1])
        self.add_rand("ref_wr_en",   [0, 1])

        # WL=16 signed integer range [-32768 .. 32767]
        self.add_rand("rx_in_re", list(range(-32768, 32768)))
        self.add_rand("rx_in_im", list(range(-32768, 32768)))
        # ref_wr_re / ref_wr_im intentionally NOT in add_rand:
        # sequences set them directly to control exact RAM contents.

    # =========================================================
    # STRING FORMATTING
    # =========================================================
    def convert2string(self):
        """Full transaction: inputs + outputs."""
        return (
            f"{self.get_name()} "
            f"IN[ rst_n={self.rst_n}, "
            f"rx_vld={self.rx_valid_in}, rx_re={self.rx_in_re}, rx_im={self.rx_in_im}, "
            f"ref_wr={self.ref_wr_en}, ref_re={self.ref_wr_re}, ref_im={self.ref_wr_im} ] -> "
            f"OUT[ ofdm_vld={self.ofdm_valid_out}, "
            f"ofdm_re={self.ofdm_out_re}, ofdm_im={self.ofdm_out_im} | "
            f"radar_vld={self.radar_valid_out}, "
            f"radar_re={self.radar_out_re}, radar_im={self.radar_out_im} ]"
        )

    def convert2string_stimulus(self):
        """Input-only: for driver / sequencer logging."""
        return (
            f"{self.get_name()} "
            f"STIMULUS: rst_n={self.rst_n}, "
            f"rx_vld={self.rx_valid_in}, rx_re={self.rx_in_re}, rx_im={self.rx_in_im}, "
            f"ref_wr={self.ref_wr_en}, ref_re={self.ref_wr_re}, ref_im={self.ref_wr_im}"
        )

    def convert2string_output(self):
        """Output-only: for monitor / scoreboard logging."""
        return (
            f"{self.get_name()} "
            f"OUTPUT: "
            f"ofdm_vld={self.ofdm_valid_out}, "
            f"ofdm_re={self.ofdm_out_re}, ofdm_im={self.ofdm_out_im} | "
            f"radar_vld={self.radar_valid_out}, "
            f"radar_re={self.radar_out_re}, radar_im={self.radar_out_im}"
        )