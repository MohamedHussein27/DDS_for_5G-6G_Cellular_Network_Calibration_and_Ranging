"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: mixed_ofdm_radar_small_test.py

    Description:
        TC-RX-005-SMALL: Micro-Signal Verification.
        Verifies the radar matched filter output using heavily constrained 
        amplitudes to definitively rule out internal FFT overflow.

    STIMULUS:
        Combined OFDM (bins 0–2047) + chirp (bins 2048–4095).
        Amplitudes artificially scaled down to a maximum integer value of 10.

    EXPECTED / PASS:
        radar_out successfully processes the frame without hitting the 
        16-bit physical ceiling. Note: Quantization noise ratio will be higher.

    Feature: F-RX-01-SMALL
"""

import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from pyuvm import *
import pyuvm

from base_test                  import base_test
from top_seq_item               import rx_item
from mixed_ofdm_radar_small_seq import mixed_ofdm_radar_small_seq


@pyuvm.test()
class mixed_ofdm_radar_small_test(base_test):
    """
    TC-RX-005-SMALL: Micro-Signal Verification.
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

        self.seq_tc_rx_005_small = mixed_ofdm_radar_small_seq.create("seq_tc_rx_005_small")
        self.seq_tc_rx_005_small.num_frames = 1
        self.seq_tc_rx_005_small.seed       = 0

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(
            f"================ Start of {self.get_type_name()} ================"
        )
        await self.generate_clock()

        self.logger.info(
            "--- Executing TC-RX-005-SMALL: Micro-Signal No-Overflow Check ---"
        )
        await self.seq_tc_rx_005_small.start(self.env.top_agt.sqr)

        # Wait for the final bit-reversal ping-pong buffer to drain
        self.logger.info("Waiting for final bit-reversal buffer to drain...")
        await ClockCycles(self.dut.clk, 2500)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT]   : {self.get_type_name()}")
        self.logger.info(" [TC ID]         : TC-RX-005-SMALL")
        self.logger.info(
            " [STIMULUS]      : Combined OFDM + chirp. Amplitudes constrained to +/- 10."
        )
        self.logger.info(
            " [PASS CRITERIA] : hardware successfully avoids internal overflow."
        )
        self.logger.info(" [FEATURE]       : F-RX-01-SMALL")
        self.logger.info("---------------------------------------------------------")