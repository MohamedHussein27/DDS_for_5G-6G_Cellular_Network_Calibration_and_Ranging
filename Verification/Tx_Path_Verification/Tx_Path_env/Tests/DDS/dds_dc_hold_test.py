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
# Test : DC Hold Test (For verifying correct handling of DC input conditions)  
# ---------------------------------------------------------                
@pyuvm.test()
class dds_dc_test(dds_base_test):
    def build_phase(self):
        super().build_phase()
        # Instantiate the DC hold sequence for this test
        self.dc_hold_seq = dds_dc_hold_seq.create("dc_hold_seq")

    async def run_phase(self):
        # Raise objection to prevent the simulation from ending prematurely
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
        
        # Start the clock and execute the hardware reset routine
        await self.generate_clock()
        await self.run_initial_setup()
        
        self.logger.info(f"Starting sequence: {self.dc_hold_seq.get_name()}") 
        
        # Execute the DC hold sequence on the DDS agent's sequencer
        await self.dc_hold_seq.start(self.env.dds_agt.sqr)
        
        self.logger.info(f"finished sequence: {self.dc_hold_seq.get_name()}")
        
        # Drop objection to allow the test to finish
        self.drop_objection() 

    def report_phase(self):
        # Print a final summary report at the end of the simulation
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : DDS output correctness")
        self.logger.info("---------------------------------------------------------")             