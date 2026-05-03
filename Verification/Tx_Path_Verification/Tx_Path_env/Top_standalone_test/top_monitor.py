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
            seq_item.rst_n          = int(self.dut_mon.rst_n.value)  
            seq_item.dds_enable     = int(self.dut_mon.dds_enable.value)
            seq_item.FTW_start      = int(self.dut_mon.FTW_start.value)
            seq_item.cycles         = int(self.dut_mon.cycles.value)
            seq_item.FTW_step       = int(self.dut_mon.FTW_step.value)


            #seq_item.ofdm_in_real   = self.dut_mon.ofdm_in_re.value
            #seq_item.ofdm_in_imag   = self.dut_mon.ofdm_in_im.value

            # capturing the IFFT outputs 
            if self.dut_mon.tx_valid.value == 1:
                seq_item.tx_valid       = int(self.dut_mon.tx_valid.value)
                seq_item.tx_out_real    = int(self.dut_mon.tx_out_re.value)
                seq_item.tx_out_imag    = int(self.dut_mon.tx_out_im.value)
            else:
                seq_item.tx_valid       = 0
                seq_item.tx_out_real    = 0
                seq_item.tx_out_imag    = 0
            
            # monitoring the top
            self.logger.info(f"Monitoring top: rst_n={seq_item.rst_n}, dds_enable={seq_item.dds_enable}, FTW_start={seq_item.FTW_start}, cycles={int(seq_item.cycles)}, FTW_step={seq_item.FTW_step}, tx_valid={seq_item.tx_valid}, tx_out_real={seq_item.tx_out_real}, tx_out_imag={seq_item.tx_out_imag}")
            self.mon_ap.write(seq_item)
            