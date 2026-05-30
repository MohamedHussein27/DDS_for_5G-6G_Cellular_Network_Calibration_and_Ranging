"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: null_input_test.py

    Description:
        TC-RX-001: Null Input — All Zeros.
        Verifies that zero input with no reference produces zero output
        on both ofdm_out and radar_out.

    STIMULUS:
        rx_in_re = 0, rx_in_im = 0 for all 4096 samples.
        rx_valid_in = 1. No reference RAM write.

    EXPECTED / PASS:
        All 2048 ofdm_out samples = 0.
        All 2048 radar_out samples = 0.
        Both valid counts correct.

    Feature: F-RX-05
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from pyuvm import *
import pyuvm

from base_test       import base_test
from top_seq_item    import rx_item
from null_input_seq  import null_input_seq


@pyuvm.test()
class null_input_test(base_test):
    """
    TC-RX-001: Null Input — All Zeros.
    All ofdm_out and radar_out samples must be zero.
    """

    def build_phase(self):
        super().build_phase()

        ConfigDB().set(self, "env.top_agt",  "is_active",
                       uvm_active_passive_enum.UVM_ACTIVE)
        ConfigDB().set(self, "env.fft_agt",  "is_active",
                       uvm_active_passive_enum.UVM_PASSIVE)
        ConfigDB().set(self, "env.ifft_agt", "is_active",
                       uvm_active_passive_enum.UVM_PASSIVE)

        ConfigDB().set(self, "env", "VERIF_MODE", "TOP")

        self.seq_tc_rx_001 = null_input_seq.create("seq_tc_rx_001")
        self.seq_tc_rx_001.num_frames = 1

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )
        await self.generate_clock()

        self.logger.info(
            "--- Executing TC-RX-001: Null Input — All Zeros ---"
        )
        await self.seq_tc_rx_001.start(self.env.top_agt.sqr)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(" [TC ID]         : TC-RX-001")
        self.logger.info(
            " [STIMULUS]      : rx_in = 0+0j for all 4096 samples. "
            "No reference RAM write."
        )
        self.logger.info(
            " [PASS CRITERIA] : "
            "All 2048 ofdm_out = 0. All 2048 radar_out = 0. "
            "Both valid counts correct."
        )
        self.logger.info(" [FEATURE]       : F-RX-05")
        self.logger.info("---------------------------------------------------------")
