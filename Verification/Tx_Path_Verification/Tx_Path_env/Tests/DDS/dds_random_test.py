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
from dds_random_sequence import * 
# ---------------------------------------------------------
# Test 2: The random test 
# ---------------------------------------------------------
@pyuvm.test()
class dds_random(dds_base_test):
    def build_phase(self):
        super().build_phase()
        self.random_seq = dds_random_seq("random_seq")
        
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        await self.generate_clock()
        await self.run_initial_setup()
        self.logger.info(f"Starting sequence: {self.random_seq.get_name()}")
        await self.random_seq.start(self.env.dds_agt.sqr)
        self.logger.info(f"Finished sequence: {self.random_seq.get_name()}")
        self.drop_objection()