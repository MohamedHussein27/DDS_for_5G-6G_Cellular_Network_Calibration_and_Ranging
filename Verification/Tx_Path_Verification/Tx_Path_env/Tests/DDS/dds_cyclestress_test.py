"""
"""
import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# Import the base test and the specific sequences
from dds_base_test import *
from dds_cyclestress_sequence import * 
# ---------------------------------------------------------
# Test 4: cyclestress test (For fast debugging of cycle counter overflow handling)
# ---------------------------------------------------------
@pyuvm.test()
class dds_cyclestress_test(dds_base_test):
    def build_phase(self):
        super().build_phase()
        self.cyclestress_seq = dds_cyclestress_seq("cyclestress_seq")
        
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.cyclestress_seq.get_name()}")
        # Runs ONLY the chirp sequence
        await self.cyclestress_seq.start(self.env.dds_agt.sqr)
        self.logger.info(f"finished sequence: {self.cyclestress_seq.get_name()}")
        self.drop_objection() 
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")    