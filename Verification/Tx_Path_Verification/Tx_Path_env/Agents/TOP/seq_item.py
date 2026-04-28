"""

"""

import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

# Multiple inheritance: UVM item features + CRV randomization features
class top_item(uvm_sequence_item, Randomized):

    def __init__(self, name="top_item"):
        # You must initialize both parent classes
        uvm_sequence_item.__init__(self, name)
        Randomized.__init__(self)

        # =========================================================
        # Input variables (Randomized)
        # =========================================================
        self.rst_n = 1
        self.dds_enable = 0
        self.FTW_start = 0
        self.cycles = 0
        self.FTW_step = 0
        self.ofdm_rd_en = 0
        self.ofdm_in_real = 0
        self.ofdm_in_imag = 0

        # =========================================================
        # Output variables (Not randomized, captured by monitor)
        # =========================================================
        self.tx_valid = 0
        self.tx_out_real = 0
        self.tx_out_imag = 0

        # =========================================================
        # CRV SETUP: Registering variables with the solver
        # =========================================================
        
        self.add_rand("rst_n", [0, 1])
        self.add_rand("dds_enable", [0, 1])
        self.FTW_start = random.randint(0, (1 << 32) - 1)
        self.cycles = random.randint(0, (1 << 13) - 1)
        self.FTW_step = random.randint(0, (1 << 32) - 1)
        
        # WL=16 Signed Integer Range: [-32768 to 32767]
        # self.add_rand("ofdm_rd_en", [1, 0])
        self.add_rand("ofdm_in_real", list(range(-32768, 32768)))  
        self.add_rand("ofdm_in_imag", list(range(-32768, 32768)))  

    # =========================================================
    # STRING FORMATTING
    # =========================================================
    def convert2string(self):
        # Full transaction printout (Inputs -> Outputs)
        return (f"{self.get_name()} "
                f"IN[ rst_n={self.rst_n}, dds_en={self.dds_enable}, "
                f"FTW_start={self.FTW_start}, cycles={self.cycles}, FTW_step={self.FTW_step}, "
                f"ofdm_rd_en={self.ofdm_rd_en}, ofdm_re={self.ofdm_in_real}, ofdm_im={self.ofdm_in_imag} ] -> "
                f"OUT[ tx_valid={self.tx_valid}, tx_re={self.tx_out_real}, tx_im={self.tx_out_imag} ]")

    def convert2string_stimulus(self):
        # Input-only printout for the driver/sequencer logging
        return (f"{self.get_name()} "
                f"STIMULUS: rst_n={self.rst_n}, dds_en={self.dds_enable}, "
                f"FTW_start={self.FTW_start}, cycles={self.cycles}, FTW_step={self.FTW_step}, "
                f"ofdm_rd_en={self.ofdm_rd_en}, ofdm_re={self.ofdm_in_real}, ofdm_im={self.ofdm_in_imag}")