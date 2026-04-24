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
from ifft_item import ifft_item

class ifft_monitor(uvm_monitor):
    def build_phase(self):
        self.dut = ConfigDB().get(self, "", "DUT")
        # Port used to broadcast observed traffic
        self.mon_port = uvm_analysis_port("mon_port", self)

    async def run_phase(self):
        while True:
            rsp_seq_item = ifft_item("rsp_seq_item")

            await RisingEdge(self.dut.clk)
            # Wait until all signals have settled for the current clock cycle
            await ReadOnly()
            
            # If the IFFT produces a valid output sample
            if self.dut.valid_out.value == 1:
                # Create a new transaction item
                rsp_seq_item = ifft_item("rsp_seq_item")
                
                # Capture the physical pin values
                rsp_seq_item.valid_out     = 1
                rsp_seq_item.data_real_out = self.dut.out_real.value.signed_integer
                rsp_seq_item.data_imag_out = self.dut.out_imag.value.signed_integer

                # capture inputs as well
                rsp_seq_item.valid_in      = self.dut.valid_in.value
                rsp_seq_item.data_real_in  = self.dut.in_real.value.signed_integer
                rsp_seq_item.data_imag_in  = self.dut.in_imag.value.signed_integer
                
                # Broadcast the captured item to the agent
                self.mon_port.write(rsp_seq_item)