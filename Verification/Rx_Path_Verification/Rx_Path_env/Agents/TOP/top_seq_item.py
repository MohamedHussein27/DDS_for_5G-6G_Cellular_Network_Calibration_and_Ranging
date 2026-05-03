"""
    Module: rx_seq_item.py

    Description:
        The sequence item defines the transaction-level packet for the RX_TOP block.
        It utilizes multiple inheritance to support both UVM reporting and native 
        Python CRV (Constrained Random Verification). It encapsulates the complex 
        datapath signals conforming to the WL=16 fixed-point format for both the
        channel data path and the reference RAM path.
"""

import cocotb 
from cocotb_coverage.crv import *
from pyuvm import * 
import random

# Multiple inheritance: UVM item features + CRV randomization features
class rx_item(uvm_sequence_item, Randomized):

    def __init__(self, name="rx_item"):
        # You must initialize both parent classes
        uvm_sequence_item.__init__(self, name)
        Randomized.__init__(self)

        # =========================================================
        # Input variables (Randomized)
        # =========================================================
        self.rst_n = 1
        
        # Channel Input Path
        self.rx_valid_in = 0
        self.rx_in_re = 0
        self.rx_in_im = 0
        
        # Reference RAM Input Path
        self.ref_wr_en = 0
        self.ref_wr_re = 0
        self.ref_wr_im = 0

        # =========================================================
        # Output variables (Not randomized, captured by monitor)
        # =========================================================
        # OFDM Output Path
        self.ofdm_valid_out = 0
        self.ofdm_out_re = 0
        self.ofdm_out_im = 0
        
        # Radar Output Path
        self.radar_valid_out = 0
        self.radar_out_re = 0
        self.radar_out_im = 0
        
        

        # =========================================================
        # CRV SETUP: Registering variables with the solver
        # =========================================================
        
        self.add_rand("rst_n", [0, 1])
        self.add_rand("rx_valid_in", [0, 1])
        self.add_rand("ref_wr_en", [0, 1])
        
        # WL=16 Signed Integer Range: [-32768 to 32767]
        self.add_rand("rx_in_re", list(range(-32768, 32768)))  
        self.add_rand("rx_in_im", list(range(-32768, 32768)))  
        #self.add_rand("ref_wr_re", list(range(-32768, 32768)))  
        #self.add_rand("ref_wr_im", list(range(-32768, 32768)))  

    # =========================================================
    # STRING FORMATTING
    # =========================================================
    def convert2string(self):
        # Full transaction printout (Inputs -> Outputs)
        return (f"{self.get_name()} "
                f"IN[ rst_n={self.rst_n}, rx_vld={self.rx_valid_in}, rx_re={self.rx_in_re}, rx_im={self.rx_in_im}, "
                f"ref_wr={self.ref_wr_en}, ref_re={self.ref_wr_re}, ref_im={self.ref_wr_im} ] -> "
                f"OUT[ ofdm_vld={self.ofdm_valid_out}, ofdm_re={self.ofdm_out_re}, ofdm_im={self.ofdm_out_im} | "
                f"radar_vld={self.radar_valid_out}, radar_re={self.radar_out_re}, radar_im={self.radar_out_im} ]")

    def convert2string_stimulus(self):
        # Input-only printout for the driver/sequencer logging
        return (f"{self.get_name()} "
                f"STIMULUS: rst_n={self.rst_n}, rx_vld={self.rx_valid_in}, rx_re={self.rx_in_re}, rx_im={self.rx_in_im}, "
                f"ref_wr={self.ref_wr_en}, ref_re={self.ref_wr_re}, ref_im={self.ref_wr_im}")