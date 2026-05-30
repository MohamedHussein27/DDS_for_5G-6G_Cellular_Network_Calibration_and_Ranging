"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: mixed_ofdm_radar_test.py

    Description:
        TC-RX-005: Mixed OFDM+Radar — OFDM Extraction.
        Verifies that ofdm_out and radar_out independently match their
        respective golden references when a combined OFDM+chirp signal
        is applied with a loaded chirp reference in the RAM.

    STIMULUS:
        Combined OFDM (bins 0–2047) + chirp (bins 2048–4095).
        Valid chirp reference loaded into reference RAM before frame.

    EXPECTED / PASS:
        ofdm_out matches golden OFDM reference.
        radar_out matches golden radar reference.
        No cross-path contamination.

    Feature: F-RX-01
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from pyuvm import *
import pyuvm

from base_test             import base_test
from top_seq_item          import rx_item
from mixed_ofdm_radar_seq  import mixed_ofdm_radar_seq


@pyuvm.test()
class mixed_ofdm_radar_test(base_test):
    """
    TC-RX-005: Mixed OFDM+Radar — OFDM Extraction.
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

        self.seq_tc_rx_005 = mixed_ofdm_radar_seq.create("seq_tc_rx_005")
        self.seq_tc_rx_005.num_frames = 1
        self.seq_tc_rx_005.seed       = 42

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )
        await self.generate_clock()

        self.logger.info(
            "--- Executing TC-RX-005: Mixed OFDM+Radar — OFDM Extraction ---"
        )
        await self.seq_tc_rx_005.start(self.env.top_agt.sqr)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(" [TC ID]         : TC-RX-005")
        self.logger.info(
            " [STIMULUS]      : Combined OFDM (bins 0–2047) + chirp "
            "(bins 2048–4095). Valid chirp reference loaded into RAM."
        )
        self.logger.info(
            " [PASS CRITERIA] : "
            "ofdm_out matches golden OFDM reference. "
            "radar_out matches golden radar reference. "
            "No cross-path contamination."
        )
        self.logger.info(" [FEATURE]       : F-RX-01")
        self.logger.info("---------------------------------------------------------")
