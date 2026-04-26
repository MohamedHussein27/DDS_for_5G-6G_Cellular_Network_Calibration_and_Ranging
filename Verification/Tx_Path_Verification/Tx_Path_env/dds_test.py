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
from base_test import base_test
from dds_sequences import * ############## need to be modified 

@pyuvm.test()
class dds_test(base_test):
    def build_phase(self):
        super().build_phase()
        ConfigDB.set(self, "*", "VERIF_MODE", "DDS")  # Set the verification mode to TOP for this test

        # write sequences here >>>>>>>

    # run phase
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        # 1. Start the clock 
        await self.generate_clock()

        # 2. Call the universal setup from base_test (Reset)
        await self.run_initial_setup()

        # ifft sequences

        # 4. Drop the objection to end the simulation
        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")