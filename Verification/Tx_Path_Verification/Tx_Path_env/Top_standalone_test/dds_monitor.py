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

class dds_monitor(uvm_monitor):
    def build_phase(self):
        self.mon_ap = uvm_analysis_port("mon_ap", self)
        self.dut_mon = ConfigDB().get(self,"","DUT")
        
    async def run_phase(self):
        sample_index = 1 # Start counter at 1 for a new chirp
        #await RisingEdge(self.dut_mon.clk)
        while True:
            await RisingEdge(self.dut_mon.clk)
            await ReadOnly()
            
            # 1. Create a new sequence item every clock cycle
            seq_item = dds_seq_item.create("seq_item")
            
            # 2. Safely capture core control signals
            try:
                seq_item.rst_n = int(self.dut_mon.rst_n.value)
            except ValueError:
                seq_item.rst_n = 0 # Assume reset if X/Z

            try:
                seq_item.enable = int(self.dut_mon.enable.value)
                seq_item.valid_out = int(self.dut_mon.valid_out.value)
            except ValueError:
                seq_item.enable = 0
                seq_item.valid_out = 0

            # 3. Capture Amplitude (handles two's complement conversion if valid)
            try:
                val = int(self.dut_mon.final_amplitude.value)
                if val > 127:
                    val -= 256
                seq_item.final_amplitude = val
            except ValueError:
                seq_item.final_amplitude = 0xFF if seq_item.rst_n == 0 else 0

            # 4. Handle state-specific logic (Reset vs Active)
            if seq_item.rst_n == 0:
                # During reset, reset the chirp tracker
                sample_index = 1
                seq_item.sample_index = 0 
            else:
                # Out of reset, safely capture streaming parameters
                try:
                    seq_item.FTW_start = int(self.dut_mon.FTW_start.value)
                    seq_item.FTW_step = int(self.dut_mon.FTW_step.value)
                    seq_item.cycles = int(self.dut_mon.cycles.value)
                except ValueError:
                    # Default to 0 if bus is X/Z during idle
                    seq_item.FTW_start = 0
                    seq_item.FTW_step = 0
                    seq_item.cycles = 0

                # Tag the sequence item with the current chirp index
                seq_item.sample_index = sample_index

                # Increment index ONLY if data is valid; otherwise reset for the next chirp
                if seq_item.valid_out == 1:
                    sample_index += 1
                else:
                    sample_index = 1
            self.logger.debug(f"Monitoring dds: rst_n={seq_item.rst_n}, enable={seq_item.enable}, FTW_start={seq_item.FTW_start}, FTW_step={seq_item.FTW_step}," 
                             f"cycles={seq_item.cycles}, valid_out={seq_item.valid_out}, "
                             f"final_amplitude={seq_item.final_amplitude}, sample_index={getattr(seq_item, 'sample_index', 'N/A')}"   )
            
            # 5. ALWAYS write to the analysis port
            self.mon_ap.write(seq_item)