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
from dds_sequences import * 
# ---------------------------------------------------------
# Test 5: Singletone only test (For fast debugging of single tone generation and FTW_step handling)
# ---------------------------------------------------------        
@pyuvm.test()
class dds_singletone_only_test(dds_base_test):
    def build_phase(self):
        super().build_phase()
        self.singletone_seq = dds_singletone_seq("singletone_seq")
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.singletone_seq.get_name()}")
        # Runs ONLY the chirp sequence
        await self.singletone_seq.start(self.env.dds_agt.sqr)
        self.logger.info(f"finished sequence: {self.singletone_seq.get_name()}")
        self.drop_objection() 
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")    