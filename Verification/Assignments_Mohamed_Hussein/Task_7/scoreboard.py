import cocotb 
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import * 
import pyuvm
import logging

class Scoreboard(uvm_scoreboard):

    def __init__(self, name, parent):
        super().__init__(name, parent)
    
    def build_phase(self):
        # Create the analysis export and the TLM FIFO
        self.sb_export = uvm_analysis_export("sb_export", self)
        self.sb_fifo = uvm_tlm_analysis_fifo("sb_fifo", self)

        # connecting here instead
        self.sb_export = self.sb_fifo.analysis_export

        # Initialize tracking counters
        self.correct_count_out = 0
        self.error_count_out = 0
        self.correct_count_c = 0
        self.error_count_c = 0

    """def connect_phase(self):
        self.sb_export.connect(self.sb_fifo.analysis_export)"""

    async def run_phase(self):
        while True:   
            seq_item_sb = await self.sb_fifo.get()

            
            out_ref, c_ref = self.ref_model(seq_item_sb)

            
            if seq_item_sb.out != out_ref:
                self.error_count_out += 1
                self.logger.error(f"Error in OUT at count {self.error_count_out} | "
                                  f"DUT OUT: {seq_item_sb.out} | REF OUT: {out_ref}")
            else:
                self.correct_count_out += 1
                
                self.logger.debug(f"Correct OUT: {seq_item_sb.out}") 

            
            if seq_item_sb.c != c_ref:
                self.error_count_c += 1
                self.logger.error(f"Error in C at count {self.error_count_c} | "
                                  f"DUT C: {seq_item_sb.c} | REF C: {c_ref}")
            else:
                self.correct_count_c += 1
                self.logger.debug(f"Correct C: {seq_item_sb.c}")
    
    
    def ref_model(self, item):

        if item.rst_n == 0:
            return 0, 0

        out_ref = 0
        c_ref = 0

        if item.op == 0:    # Addition
            res = item.a + item.b
            out_ref = res & 0xF     # Mask to 4 bits (equivalent to 4'b)
            c_ref = (res >> 4) & 1  # Shift right by 4 to grab the 5th bit (carry)
            
        elif item.op == 1:  # XOR
            out_ref = item.a ^ item.b
            c_ref = 0
            
        elif item.op == 2:  # AND
            out_ref = item.a & item.b
            c_ref = 0
            
        elif item.op == 3:  # OR
            out_ref = item.a | item.b
            c_ref = 0

        return out_ref, c_ref

    def report_phase(self):
        # Print the final statistics at the end of the test
        self.logger.info("================ SCOREBOARD REPORT ================")
        self.logger.info(f" Total successful OUT transactions: {self.correct_count_out}")
        self.logger.info(f" Total failed OUT transactions:     {self.error_count_out}")
        self.logger.info(f" Total successful C transactions:   {self.correct_count_c}")
        self.logger.info(f" Total failed C transactions:       {self.error_count_c}")
        self.logger.info("===================================================")