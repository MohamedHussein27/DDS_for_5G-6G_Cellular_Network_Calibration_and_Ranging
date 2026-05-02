"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
    
    Module: base_test.py

    Description:
        The base test serves as the foundational pyuvm test class for the TX wrapper 
        verification environment. It establishes the testbench hierarchy, instantiates 
        the environment components, generates the 491.52 MHz system clock, and defines 
        the standard execution flow (reset followed by the main sequence) that all 
        extended tests will inherit.
"""


import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging


from environment import *
from sequences import *


class base_test(uvm_test):

    def __init__(self, name, parent):
        super().__init__(name, parent)


    def build_phase(self):
        
        # 1. Environment 
        self.env = Environment.create("env", self)

        # dut handle
        self.dut = cocotb.top

        # put dut handle to be seen by all components
        ConfigDB().set(self, "*", "DUT", self.dut)

        
        self.reset_seq = reset_sequence.create("reset_seq")
        self.seq       = base_sequence.create("seq")

    # generate the clock
    async def generate_clock(self):
        self.clk = Clock(self.dut.clk, 2034.5, units="ps") # 491.52 MHz
        await cocotb.start(self.clk.start())


    async def run_phase(self):
        self.raise_objection()

        # 1. Start the clock using the task prepared in build_phase
        await self.generate_clock()

        # 2. Run Reset
        self.logger.info(f"********** Start of {self.get_type_name()} **********")
        self.logger.info("Starting Reset Sequence...")
        await self.reset_seq.start(self.env.agt.sqr)

        # 3. Run Main Sequence 
        self.logger.info(f"Starting sequence: {self.seq.get_type_name()}")
        await self.seq.start(self.env.agt.sqr)
        self.logger.info(f"Finished sequence: {self.seq.get_type_name()}")

        self.drop_objection()
        
    def final_phase(self):
        # Prints the factory overrides to the transcript to prove the swap worked
        uvm_factory().print(0)
        self.logger.info(f"********** End of {self.get_type_name()} **********")