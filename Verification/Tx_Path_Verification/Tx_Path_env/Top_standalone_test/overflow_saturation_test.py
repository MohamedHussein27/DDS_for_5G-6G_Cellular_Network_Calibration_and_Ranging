"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: overflow_saturation_test.py

    Description:
        This test verifies the TX datapath under TC-007: Overflow / Saturation Behaviour.
        It drives worst-case maximum-amplitude stimulus (max M_dds + all OFDM corner
        symbols simultaneously).
        If any sample hits ±32767, the next sample must NOT jump to the opposite extreme
        (no wrap-around). Saturate-don't-wrap: verify sign does not invert at peak amplitude.
        Feature: F-02
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging

# Import the base test and the specific sequences
from base_test import base_test
from top_seq_item import *
from overflow_saturation_seq import overflow_saturation_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_007_overflow_saturation_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_overflow_saturation = overflow_saturation_seq.create("seq_overflow_saturation")
        self.seq_reset_before        = reset_before_frame_seq.create("seq_reset_before")

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the Overflow / Saturation Sequence
        self.logger.info("--- Executing TC-007: Overflow / Saturation Behaviour ---")
        await self.seq_overflow_saturation.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Worst-case max-amplitude, all OFDM corner symbols simultaneously")
        self.logger.info(" [PASS CRITERIA]: No wrap-around on any sample. Sign must NOT invert at peak amplitude.")
        self.logger.info(" [FEATURE]     : F-02 — Saturate-don't-wrap: verify sign does not invert at peak amplitude")
        self.logger.info("---------------------------------------------------------")
