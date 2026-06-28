"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: output_word_format_test.py

    Description:
        This test verifies the TX datapath under TC-006: Output Word Format — Normal Operation.
        It drives a dual-band nominal stimulus and inspects all 4096 tx_out samples for
        correct signed 16-bit fixed-point format (WL=16, FL=5).
        All samples must be within [-32768, +32767].
        No sample equals -32768 unless saturation is explicitly expected.
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
from output_word_format_seq import output_word_format_seq
from reset_sequence import reset_before_frame_seq

@pyuvm.test()
class tc_006_output_word_format_test(base_test):
    def build_phase(self):
        super().build_phase()

        # Instantiate only the sequences needed for this test
        self.seq_output_word_format = output_word_format_seq.create("seq_output_word_format")
        self.seq_reset_before       = reset_before_frame_seq.create("seq_reset_before")

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock
        await self.generate_clock()

        self.logger.info("--- Executing TC-009: Reset Before Frame Start ---")
        await self.seq_reset_before.start(self.env.top_agt.sqr)

        # 2. Execute the Output Word Format Sequence
        self.logger.info("--- Executing TC-006: Output Word Format — Normal Operation ---")
        await self.seq_output_word_format.start(self.env.top_agt.sqr)

        # 3. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Signed 16-bit fixed-point format (WL=16, FL=5)")
        self.logger.info(" [PASS CRITERIA]: All samples in [-32768, +32767]. No -32768 unless saturation expected.")
        self.logger.info(" [FEATURE]     : F-02 — Format check: signed 16-bit two's complement, FL=5 fractional bits")
        self.logger.info("---------------------------------------------------------")
