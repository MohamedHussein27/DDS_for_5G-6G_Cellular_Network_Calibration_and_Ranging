"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: monitor.py
"""
import pyuvm
from pyuvm import *
from dds_seq_item import *
import cocotb
from cocotb.triggers import RisingEdge, ReadOnly
from dds_seq_item import dds_seq_item

class dds_monitor(uvm_monitor):
    def build_phase(self):
        self.mon_ap = uvm_analysis_port("mon_ap", self)
        self.dut_mon = ConfigDB().get(self,"","DUT")
        
    async def run_phase(self):
        sample_index = 1 # Start counter at 1 for a new chirp
        await RisingEdge(self.dut_mon.clk)
        while True:
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            # 1. Continuous Reset Checking
            if self.dut_mon.rst_n.value == 0:
                seq_item = dds_seq_item.create("seq_item")
                seq_item.rst_n = 0
                seq_item.enable = int(self.dut_mon.enable.value)
                seq_item.valid_out = int(self.dut_mon.valid_out.value)
                self.logger.info(f"Monitoring dds: rst_n={seq_item.rst_n}, enable={seq_item.enable}, FTW_start={seq_item.FTW_start}, FTW_step={seq_item.FTW_step}, cycles={seq_item.cycles}, valid_out={seq_item.valid_out}, final_amplitude={seq_item.final_amplitude}, sample_index={getattr(seq_item, 'sample_index', 'N/A')}"   )
                try:
                    seq_item.final_amplitude = int(self.dut_mon.final_amplitude.value)
                except ValueError:
                    seq_item.final_amplitude = 0xFF 
                
                self.mon_ap.write(seq_item)
                sample_index = 1 # Reset the tracker while in reset
                
                while self.dut_mon.rst_n.value == 0:
                    await RisingEdge(self.dut_mon.clk)
                    await ReadOnly()
                continue
            
            # 2. Continuous Data Streaming
            if self.dut_mon.valid_out.value == 1:
                seq_item = dds_seq_item.create("seq_item")
                self.logger.info(f"Monitoring dds: rst_n={seq_item.rst_n}, enable={seq_item.enable}, FTW_start={seq_item.FTW_start}, FTW_step={seq_item.FTW_step}, cycles={seq_item.cycles}, valid_out={seq_item.valid_out}, final_amplitude={seq_item.final_amplitude}, sample_index={getattr(seq_item, 'sample_index', 'N/A')}"   )
                seq_item.rst_n = 1
                seq_item.FTW_start = int(self.dut_mon.FTW_start.value)
                seq_item.FTW_step = int(self.dut_mon.FTW_step.value)
                seq_item.cycles = int(self.dut_mon.cycles.value)
                seq_item.enable = int(self.dut_mon.enable.value)
                seq_item.valid_out = 1
                seq_item.sample_index = sample_index # Tag the specific cycle
                
                val = int(self.dut_mon.final_amplitude.value)
                if val > 127:
                    val -= 256
                seq_item.final_amplitude = val
                
                # Send the single cycle instantly!
                self.mon_ap.write(seq_item)
                sample_index += 1
            else:
                # If valid_out drops, the hardware chirp ended. Reset for the next one.
                sample_index = 1

            