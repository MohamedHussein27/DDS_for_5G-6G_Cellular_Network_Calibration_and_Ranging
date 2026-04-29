"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: fft_monitor.py

    Description:
        The monitor is a passive component that observes the FFT outputs. 
        It constantly watches the valid_out signal. Whenever valid_out is asserted 
        (high), it captures the out_real and out_imag values, packages them into an 
        fft_item, and broadcasts the item through its analysis port to the 
        scoreboard and coverage subscribers.
"""

import cocotb
from cocotb.triggers import *
from pyuvm import *
from fft_seq_item import fft_item

class fft_monitor(uvm_monitor):
    def build_phase(self):
        self.dut = ConfigDB().get(self, "", "DUT")
        # Port used to broadcast observed traffic
        self.mon_port = uvm_analysis_port("mon_port", self)

    async def run_phase(self):
        # Wait for initial reset state to settle before entering the main loop
        await RisingEdge(self.dut.clk)
        
        while True:
            await RisingEdge(self.dut.clk)
            # Wait
            await ReadOnly()
            
            # Create a new transaction item
            rsp_seq_item = fft_item("rsp_seq_item")
            
            # ─────────────────────────────────────────────────────────
            # Safely capture control signals 
            # ─────────────────────────────────────────────────────────
            try:
                rsp_seq_item.rst_n = int(self.dut.rst_n.value)
            except ValueError:
                rsp_seq_item.rst_n = 0
                
            try:
                rsp_seq_item.valid_in = int(self.dut.valid_in.value)
            except ValueError:
                rsp_seq_item.valid_in = 0
                
            try:
                rsp_seq_item.valid_out = int(self.dut.valid_out.value)
            except ValueError:
                rsp_seq_item.valid_out = 0

            # ─────────────────────────────────────────────────────────
            # Capture Input Data (Only if valid to prevent X-propagation crashes)
            # ─────────────────────────────────────────────────────────
            if rsp_seq_item.valid_in == 1:
                rsp_seq_item.in_real = self.dut.in_real.value.signed_integer
                rsp_seq_item.in_imag = self.dut.in_imag.value.signed_integer
            else:
                rsp_seq_item.in_real = 0
                rsp_seq_item.in_imag = 0

            # ─────────────────────────────────────────────────────────
            # Capture Output Data
            # ─────────────────────────────────────────────────────────
            if rsp_seq_item.valid_out == 1:
                rsp_seq_item.out_real = self.dut.out_real.value.signed_integer
                rsp_seq_item.out_imag = self.dut.out_imag.value.signed_integer
            else:
                rsp_seq_item.out_real = 0
                rsp_seq_item.out_imag = 0

            # monitoring captured output values (uncomment for debug)
            # self.logger.info(f"Monitor captured: valid_out={rsp_seq_item.valid_out}, out_real={rsp_seq_item.out_real}, out_imag={rsp_seq_item.out_imag}")
            
            # Broadcast the captured item to the agent
            self.mon_port.write(rsp_seq_item)