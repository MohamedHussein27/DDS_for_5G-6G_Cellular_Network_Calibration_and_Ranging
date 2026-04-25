"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
    
    Module: base_test.py

    Description:
        The base test serves as the foundational pyuvm test class for the TX wrapper 
        verification environment. It establishes the testbench hierarchy, instantiates 
        the environment components, generates the 491.52 MHz system clock, and defines 
        the standard execution flow (reset only and creating the clock) that all 
        extended tests will inherit.
"""


import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm


from environment import *
from sequences import *


class base_test(uvm_test):

    def __init__(self, name, parent):
        super().__init__(name, parent)


    def build_phase(self):
        
        # Environment 
        self.env = Environment.create("env", self)

        # dut handle
        self.dut = cocotb.top

        # put dut handle to be seen by all components
        ConfigDB().set(self, "*", "DUT", self.dut)

        self.reset_seq = reset_sequence.create("reset_seq")

        # assigning the which is active and which is passive
        ConfigDB().set(self, "top_agent", "is_active", UVM_ACTIVE) # top agent
        ConfigDB().set(self, "dds_agent", "is_active", UVM_PASSIVE) # dds agent
        ConfigDB().set(self, "fft_agent", "is_active", UVM_PASSIVE) # fft agent
        ConfigDB().set(self, "ifft_agent", "is_active", UVM_PASSIVE) # ifft agent

    # generate the clock
    async def generate_clock(self):
        self.clk = Clock(self.dut.clk, 2034.5, units="ps") # 491.52 MHz
        await cocotb.start(self.clk.start())

    # reset
    async def run_initial_setup(self):
        self.logger.info("Starting Reset Sequence...")
        await self.reset_seq.start(self.env.agt.sqr)
        self.logger.info("Reset Sequence Complete.")
        
    def final_phase(self):
        self.logger.info(f"********** End of {self.get_type_name()} **********")