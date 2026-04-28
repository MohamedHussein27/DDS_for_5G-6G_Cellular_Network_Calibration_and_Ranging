"""
 
"""

import pyuvm
from pyuvm import *
from seq_item import *
import cocotb
from cocotb.clock import Clock
from cocotb.triggers import *
class top_monitor(uvm_monitor):
    def build_phase(self):
        self.mon_port = uvm_analysis_port("mon_port", self)
        self.dut_mon= ConfigDB().get(self,"","DUT")
        
    async def run_phase(self):
        while True:
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            seq_item = seq_item.create("seq_item")
            # here we will read the signals from the DUT and write them to the analysis port
            seq_item.rst_n = self.dut.rst_n.value
            seq_item.dds_enable = self.dut.dds_enable.value
            seq_item.FTW_start = self.dut.FTW_start.value
            seq_item.cycles = self.dut.cycles.value
            seq_item.FTW_step = self.dut.FTW_step.value
            seq_item.ofdm_rd_en = self.dut.ofdm_rd_en.value
            seq_item.ofdm_in_real=self.dut.ofdm_in_real.value
            seq_item.ofdm_in_imag=self.dut.ofdm_in_imag.value
            # capturing the IFFT outputs 
            seq_item.tx_valid=self.dut.tx_valid.value
            seq_item.tx_out_real=self.dut.tx_out_real.value
            seq_item.tx_out_imag=self.dut.tx_out_imag.value
            
            
            
            
            self.mon_port.write(seq_item)
            