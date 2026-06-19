"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: mixed_ofdm_radar_float_test.py

    Description:
        TC-RX-005-FLOAT: Mixed OFDM+Radar — Floating Point Verification.
"""

import cocotb
from cocotb.triggers import *
from pyuvm import *
import pyuvm

from base_test                  import base_test
from top_seq_item               import rx_item
from mixed_ofdm_radar_float_seq import mixed_ofdm_radar_float_seq


@pyuvm.test()
class mixed_ofdm_radar_float_test(base_test):
    """
    TC-RX-005-FLOAT: Mixed OFDM+Radar — Floating Point Verification.
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

        self.seq_tc_rx_005_float = mixed_ofdm_radar_float_seq.create("seq_tc_rx_005_float")
        self.seq_tc_rx_005_float.num_frames = 1
        self.seq_tc_rx_005_float.seed       = 0

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )
        await self.generate_clock()

        self.logger.info(
            "--- Executing TC-RX-005-FLOAT: Mixed OFDM+Radar (Unquantized Float) ---"
        )
        await self.seq_tc_rx_005_float.start(self.env.top_agt.sqr)

        # Wait for the final bit-reversal ping-pong buffer to drain
        self.logger.info("Waiting for final bit-reversal buffer to drain...")
        await ClockCycles(self.dut.clk, 2500)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(" [TC ID]         : TC-RX-005-FLOAT")
        self.logger.info(
            " [STIMULUS]      : Combined OFDM + chirp. Unquantized Float Path."
        )
        self.logger.info(
            " [PASS CRITERIA] : radar_out imaginary peak evaluates to exactly 0.0."
        )
        self.logger.info("---------------------------------------------------------")