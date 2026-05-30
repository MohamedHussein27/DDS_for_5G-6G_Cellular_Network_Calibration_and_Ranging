"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: rtl_hex_chirp_test.py

    Description:
        This test verifies the RX datapath under TC-RX-003:
        RTL Hex File Stimulus, Chirp Signal In.

        The reference RAM is pre-loaded from the RTL-generated hex files
        (rtl_rx_radar_re.hex / rtl_rx_radar_im.hex).
        One full 4096-sample rx_in frame is driven from the RTL-generated
        hex files (rtl_rx_in_re.hex / rtl_rx_in_im.hex).

    STIMULUS:
        Reference RAM ← rtl_rx_radar_re.hex / rtl_rx_radar_im.hex (2048 entries)
        rx_in frame   ← rtl_rx_in_re.hex   / rtl_rx_in_im.hex   (4096 samples)
        All values are 16-bit signed two's-complement (Q8.8).

    EXPECTED / PASS:
        radar_out and ofdm_out match the golden RTL outputs paired with
        the hex stimulus files.

    Feature: F-RX-03
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from pyuvm import *
import pyuvm
import logging

from base_test            import base_test
from rx_top_seq_item         import *
from rtl_hex_chirp_seq    import rtl_hex_chirp_seq


@pyuvm.test()
class rtl_hex_chirp_test(base_test):
    """
    TC-RX-003: RTL Hex File Stimulus, Chirp Signal In.
    Verifies that the RX datapath produces the expected radar_out and
    ofdm_out when driven with the RTL-generated hex file stimulus.
    """

    def build_phase(self):
        super().build_phase()

        # Override agent activity for RX TOP:
        #   rx_agt  → ACTIVE  (drives rx_valid_in + ref_wr_en)
        #   fft_agt  → PASSIVE (monitors FFT boundary)
        #   ifft_agt → PASSIVE (monitors IFFT boundary)
        ConfigDB().set(self, "env.rx_agt",  "is_active",
                       uvm_active_passive_enum.UVM_ACTIVE)
        #ConfigDB().set(self, "env.fft_agt",  "is_active",
         #              uvm_active_passive_enum.UVM_PASSIVE)
        #ConfigDB().set(self, "env.ifft_agt", "is_active",
         #              uvm_active_passive_enum.UVM_PASSIVE)

        # Set full-system verification mode
        ConfigDB().set(self, "env", "VERIF_MODE", "RX")

        # Instantiate the TC-RX-003 sequence
        self.seq_tc_rx_003 = rtl_hex_chirp_seq.create("seq_tc_rx_003")

        # Configure: run one frame (extend to stress-test by setting > 1)
        self.seq_tc_rx_003.num_frames = 1

        # Optional: override hex file paths if they are not co-located
        # with the test. Leave commented out to use the defaults.
        #
        # self.seq_tc_rx_003.rx_in_re_path  = "/path/to/rtl_rx_in_re.hex"
        # self.seq_tc_rx_003.rx_in_im_path  = "/path/to/rtl_rx_in_im.hex"
        # self.seq_tc_rx_003.ref_re_path    = "/path/to/rtl_rx_radar_re.hex"
        # self.seq_tc_rx_003.ref_im_path    = "/path/to/rtl_rx_radar_im.hex"

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )

        # 1. Start the 491.52 MHz system clock
        await self.generate_clock()

        # 2. Execute TC-RX-003
        self.logger.info(
            "--- Executing TC-RX-003: RTL Hex File Stimulus, Chirp Signal In ---"
        )
        await self.seq_tc_rx_003.start(self.env.rx_agt.sqr)

        self.drop_objection()

    def report_phase(self):
        self.logger.info(
            "---------------------------------------------------------"
        )
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(
            " [TC ID]         : TC-RX-003"
        )
        self.logger.info(
            " [STIMULUS]      : rx_in from rtl_rx_in_re/im.hex (4096 samples). "
            "Reference RAM from rtl_rx_radar_re/im.hex (2048 entries)."
        )
        self.logger.info(
            " [PASS CRITERIA] : "
            "radar_out and ofdm_out match RTL golden outputs."
        )
        self.logger.info(
            " [FEATURE]       : F-RX-03 — "
            "Full datapath verification using RTL-generated hex stimulus."
        )
        self.logger.info(
            "---------------------------------------------------------"
        )
