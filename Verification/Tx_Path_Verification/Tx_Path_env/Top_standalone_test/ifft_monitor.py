"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_monitor.py

    Description:
        The monitor is a passive component that observes the IFFT outputs. 
        It constantly watches the valid_out signal. Whenever valid_out is asserted 
        (high), it captures the out_real and out_imag values, packages them into an 
        ifft_item, and broadcasts the item through its analysis port to the 
        scoreboard and coverage subscribers.
"""

import cocotb
from cocotb.triggers import *
from pyuvm import *
from ifft_seq_item import ifft_item
from ifft_seq_item import *

def safe_int(val, default=0):
    try: return int(val)
    except ValueError: return default

def safe_signed(val, default=0):
    try: return val.signed_integer
    except ValueError: return default

class ifft_monitor(uvm_monitor):
    def build_phase(self):


        self.dut = ConfigDB().get(self, "", "DUT")
        self.mon_port = uvm_analysis_port("mon_port", self)

    async def run_phase(self):
        await RisingEdge(self.dut.clk)
        
        while True:
            await RisingEdge(self.dut.clk)
            await ReadOnly()
            
            rsp_seq_item = ifft_item("rsp_seq_item")
            
            # Safe capture to prevent X/Z propagation crashes
            rsp_seq_item.valid_out     = safe_int(self.dut.valid_out.value)
            rsp_seq_item.out_real = safe_signed(self.dut.out_real.value)
            rsp_seq_item.out_imag = safe_signed(self.dut.out_imag.value)

            rsp_seq_item.rst_n         = safe_int(self.dut.rst_n.value)
            rsp_seq_item.valid_in      = safe_int(self.dut.valid_in.value)
            rsp_seq_item.in_real  = safe_signed(self.dut.in_real.value)
            rsp_seq_item.in_imag  = safe_signed(self.dut.in_imag.value)
            
            self.logger.info(f"Monitor captured: rst_n={rsp_seq_item.rst_n}, valid_in={rsp_seq_item.valid_in},"
                              f" in_real={rsp_seq_item.in_real}, in_imag={rsp_seq_item.in_imag}, valid_out={rsp_seq_item.valid_out},"
                               f" out_real={rsp_seq_item.out_real}, out_imag={rsp_seq_item.out_imag}")
            self.mon_port.write(rsp_seq_item)