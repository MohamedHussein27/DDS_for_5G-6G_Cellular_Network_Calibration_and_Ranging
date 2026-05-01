"""
 
"""

import pyuvm
from pyuvm import *
from top_seq_item import *
import cocotb
from cocotb.clock import Clock
from cocotb.triggers import *
class top_monitor(uvm_monitor):
    def build_phase(self):
        self.mon_ap = uvm_analysis_port("mon_ap", self)
        self.dut_mon= ConfigDB().get(self,"","DUT")
        
    async def run_phase(self):
        await RisingEdge(self.dut_mon.clk)
        while True:
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            seq_item = top_item.create("seq_item")
            # here we will read the signals from the DUT and write them to the analysis port
            seq_item.rst_n = self.dut_mon.rst_n.value
            seq_item.dds_enable = self.dut_mon.dds_enable.value
            seq_item.FTW_start = self.dut_mon.FTW_start.value
            seq_item.cycles = self.dut_mon.cycles.value
            seq_item.FTW_step = self.dut_mon.FTW_step.value
            seq_item.ofdm_rd_en = self.dut_mon.ofdm_rd_en.value
            seq_item.ofdm_in_real=self.dut_mon.ofdm_in_re.value
            seq_item.ofdm_in_imag=self.dut_mon.ofdm_in_im.value
            # capturing the IFFT outputs 
            seq_item.tx_valid=self.dut_mon.tx_valid.value
            seq_item.tx_out_real=self.dut_mon.tx_out_real.value
            seq_item.tx_out_imag=self.dut_mon.tx_out_imag.value
            
            
            
            
            self.mon_ap.write(seq_item)
            