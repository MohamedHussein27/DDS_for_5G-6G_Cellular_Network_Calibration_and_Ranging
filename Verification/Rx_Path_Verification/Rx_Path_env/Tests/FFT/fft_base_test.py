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
#from sequences import *
from fft_reset_seq import FftResetSeq



class base_test(uvm_test):

    def __init__(self, name, parent):
        super().__init__(name, parent)

    def build_phase(self):
        super().build_phase()
        
        # 1. Environment instantiation
        self.env = Environment.create("env", self)

        # 2. DUT handle
        self.dut = cocotb.top

        # 3. Share DUT handle via ConfigDB
        ConfigDB().set(self, "*", "DUT", self.dut)

        # 4. Initialize Reset Sequence
        self.reset_seq = FftResetSeq.create("reset_seq")

        # 5. Configure Agent topology
        # Using the full enum path to avoid any import confusion
        ConfigDB().set(self, "fft_agent", "is_active", pyuvm.uvm_active_passive_enum.UVM_ACTIVE)

    # --- THE CLOCK FIX ---
    async def generate_clock(self):
        """
        Generates the 491.52 MHz system clock.
        Rounded to 2034ps to match the simulator's 1ps precision grid.
        """
        # Changed 2034.5 to 2035 and units="ps" to unit="ps"
        self.clk = Clock(self.dut.clk, 2034, unit="ps") 
        cocotb.start_soon(self.clk.start())

    # --- RESET ---
    async def run_initial_setup(self):
        self.logger.info("Starting Reset Sequence...")
        # Ensure the sequence starts on the FFT sequencer
        await self.reset_seq.start(self.env.fft_agt.sqr)
        self.logger.info("Reset Sequence Complete.")
        
    def final_phase(self):
        self.logger.info(f"********** End of {self.get_type_name()} **********")