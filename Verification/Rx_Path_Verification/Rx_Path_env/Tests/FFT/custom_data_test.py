"""
    Module: custom_data_test.py
"""
import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *
import pyuvm

# Component Imports
from fft_base_test import base_test 
from custom_data_seq import SeqCustomData

@pyuvm.test()
class test_fft_custom_data(base_test):
    """
    Test Case: Custom User Data Streaming
    Description: Configures the environment in FFT-only mode and runs 
    the data provided in rtl_dds_out_testing.txt.
    """
    def build_phase(self):
        # Force the environment to build the isolated FFT path
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Universal Setup
        await self.generate_clock()
        await self.run_initial_setup()
        
        self.logger.info("================================================")
        self.logger.info(" Starting Custom Data FFT Test...               ")
        self.logger.info("================================================")
        
        # 2. Start the custom data sequence
        seq = SeqCustomData("custom_seq")
        # Attach the active sequencer
        seq.my_sqr = self.env.fft_agt.sqr 
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Final Pipeline Flush
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        # 4. End Simulation
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()