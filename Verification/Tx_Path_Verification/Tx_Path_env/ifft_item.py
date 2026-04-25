"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_item.py

    Description:
        The sequence item defines the transaction-level packet for the IFFT block.
        It utilizes multiple inheritance to support both UVM reporting and native 
        Python CRV (Constrained Random Verification). It encapsulates the complex 
        datapath signals conforming to the WL=16 fixed-point format.
"""

import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# Multiple inheritance: UVM item features + CRV randomization features
class ifft_item(uvm_sequence_item, Randomized):

    def __init__(self, name="ifft_item"):
        # You must initialize both parent classes
        uvm_sequence_item.__init__(self, name)
        Randomized.__init__(self)

        # =========================================================
        # Input variables (Randomized)
        # =========================================================
        self.rst_n = 1
        self.valid_in = 0
        self.in_real = 0
        self.in_imag = 0

        # =========================================================
        # Output variables (Not randomized, captured by monitor)
        # =========================================================
        self.valid_out = 0
        self.out_real = 0
        self.out_imag = 0

        # =========================================================
        # CRV SETUP: Registering variables with the solver
        # =========================================================
        
        self.add_rand("rst_n", [0, 1])
        self.add_rand("valid_in", [0, 1])
        
        # WL=16 Signed Integer Range: [-32768 to 32767]
        self.add_rand("in_real", list(range(-32768, 32768)))  
        self.add_rand("in_imag", list(range(-32768, 32768)))  

    # =========================================================
    # STRING FORMATTING
    # =========================================================
    def convert2string(self):
        # Full transaction printout (Inputs -> Outputs)
        return (f"{self.get_name()} "
                f"IN[ rst_n={self.rst_n}, valid={self.valid_in}, real={self.in_real}, imag={self.in_imag} ] -> "
                f"OUT[ valid={self.valid_out}, real={self.out_real}, imag={self.out_imag} ]")

    def convert2string_stimulus(self):
        # Input-only printout for the driver/sequencer logging
        return (f"{self.get_name()} "
                f"STIMULUS: rst_n={self.rst_n}, valid_in={self.valid_in}, real={self.in_real}, imag={self.in_imag}")