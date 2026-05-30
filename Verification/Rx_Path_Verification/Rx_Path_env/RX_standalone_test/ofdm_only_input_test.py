"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ofdm_only_input_test.py

    Description:
        TC-RX-003: OFDM-Only Input (Radar Bins = 0).
        Verifies that an OFDM-only signal (bins 0–2047) passes through
        the OFDM path correctly and produces zero on the radar path.

    STIMULUS:
        OFDM signal in bins 0–2047. Bins 2048–4095 = 0.
        Reference RAM = zeros (never written).

    EXPECTED / PASS:
        ofdm_out matches golden model.
        radar_out = 0 on all 2048 samples.
        No cross-path contamination.

    Feature: F-RX-01
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from pyuvm import *
import pyuvm

from base_test            import base_test
from top_seq_item         import rx_item
from ofdm_only_input_seq  import ofdm_only_input_seq


@pyuvm.test()
class ofdm_only_input_test(base_test):
    """
    TC-RX-003: OFDM-Only Input (Radar Bins = 0).
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

        self.seq_tc_rx_003 = ofdm_only_input_seq.create("seq_tc_rx_003")
        self.seq_tc_rx_003.num_frames = 1
        self.seq_tc_rx_003.seed       = 42

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )
        await self.generate_clock()

        self.logger.info(
            "--- Executing TC-RX-003: OFDM-Only Input (Radar Bins = 0) ---"
        )
        await self.seq_tc_rx_003.start(self.env.top_agt.sqr)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(" [TC ID]         : TC-RX-003")
        self.logger.info(
            " [STIMULUS]      : OFDM signal in bins 0–2047. "
            "Radar bins 2048–4095 = 0. Reference RAM = zeros."
        )
        self.logger.info(
            " [PASS CRITERIA] : "
            "ofdm_out matches golden model. "
            "radar_out = 0 on all 2048 samples. No cross-path contamination."
        )
        self.logger.info(" [FEATURE]       : F-RX-01")
        self.logger.info("---------------------------------------------------------")
