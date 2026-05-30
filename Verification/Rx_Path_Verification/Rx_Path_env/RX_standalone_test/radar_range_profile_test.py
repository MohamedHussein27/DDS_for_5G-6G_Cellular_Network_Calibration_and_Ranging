"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: radar_range_profile_test.py

    Description:
        TC-RX-007: Radar Range Profile — Loaded Reference.
        Verifies bit-exact match between radar_out and the golden
        range-profile reference when a matched TX chirp reference is
        pre-loaded and a matching received chirp is applied.

    STIMULUS:
        Reference RAM ← TX chirp reference (2048 entries).
        rx_in         ← matching received chirp (4096 samples).

    EXPECTED / PASS:
        Bit-exact match on all 2048 radar_out samples (Error = 0).

    Feature: F-RX-01, F-RX-02
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from pyuvm import *
import pyuvm

from base_test                import base_test
from top_seq_item             import rx_item
from radar_range_profile_seq  import radar_range_profile_seq


@pyuvm.test()
class radar_range_profile_test(base_test):
    """
    TC-RX-007: Radar Range Profile — Loaded Reference.
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

        self.seq_tc_rx_007 = radar_range_profile_seq.create("seq_tc_rx_007")
        self.seq_tc_rx_007.num_frames = 1

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )
        await self.generate_clock()

        self.logger.info(
            "--- Executing TC-RX-007: Radar Range Profile — Loaded Reference ---"
        )
        await self.seq_tc_rx_007.start(self.env.top_agt.sqr)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(" [TC ID]         : TC-RX-007")
        self.logger.info(
            " [STIMULUS]      : TX chirp reference loaded into RAM. "
            "Matching received chirp applied as rx_in."
        )
        self.logger.info(
            " [PASS CRITERIA] : "
            "Bit-exact match on all 2048 radar_out samples (Error = 0)."
        )
        self.logger.info(" [FEATURE]       : F-RX-01, F-RX-02")
        self.logger.info("---------------------------------------------------------")
