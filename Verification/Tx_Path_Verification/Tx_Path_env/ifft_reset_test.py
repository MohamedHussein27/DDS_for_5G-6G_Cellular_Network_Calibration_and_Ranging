"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: ifft_reset_test.py
 
    Description:
        Asserts rst_n = 0 and verifies that all delay registers flush
        and valid_out stays de-asserted. No frequency-domain data is sent.
 
    Sequences
    ─────────
        reset_sequence  – drives rst_n = 0 for one item
"""
 
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging
 
from ifft_base_test  import ifft_base_test
from ifft_sequences  import reset_sequence
 
 
@pyuvm.test()
class ifft_reset_test(ifft_base_test):
 
    def build_phase(self):
        super().build_phase()
        self.seq = reset_sequence.create("seq")
 
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
 
        await self.generate_clock()
        await self.run_initial_setup()
 
        self.logger.info("Running: reset_sequence")
        await self._run(self.seq)
 
        self.drop_objection()
 
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : rst_n flushes all delay registers")
        self.logger.info("                 valid_out stays de-asserted throughout")
        self.logger.info("---------------------------------------------------------")