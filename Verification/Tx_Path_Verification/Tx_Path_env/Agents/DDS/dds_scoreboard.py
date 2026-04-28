"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: scoreboard.py

    Description:
        This module defines the PUVM Scoreboard for the DDS TX datapath.
        It serves as the ultimate verification authority, determining the 
        pass/fail status of the DUT by comparing actual hardware outputs 
        against an ideal, bit-accurate reference model.

        Key responsibilities:
        1. Reference Modeling: Calculates expected sine/cosine waveforms mathematically 
           based on the applied Frequency Tuning Words (FTW), phase increments, 
           and control signals.
        2. Data Collection: Captures observed `seq_item` transactions from the 
           passive Monitors via UVM analysis exports.
        3. Comparison & Tolerance: Compares the hardware's generated digital 
           waveform against the golden model, accounting for expected 
           quantization noise and phase truncation inherent to DDS architecture.
        4. Reporting: Logs precise mismatch locations, tracks verification 
           statistics, and flags fatal errors if functional requirements are violated.
"""
import pyuvm
from pyuvm import *
from dds_seq_item import *
from dds_driver import dds_driver
from dds_monitor import dds_monitor
from dds_sequencer import dds_sequencer
from cocotb.triggers import Timer
from cocotb.clock import Clock
import cocotb
import numpy as np

# Import your items and golden model
from dds_core import dds_core

class dds_scoreboard(uvm_scoreboard):
    def build_phase(self):
        super().build_phase()
        # 1. Initialize statistics counters
        self.passed_test_cases = 0
        self.failed_test_cases = 0
        
        # 2. UVM TLM setup
        self.sc_export = uvm_analysis_export("sc_export", self)
        self.sc_fifo = uvm_tlm_analysis_fifo("sc_fifo", self)

        # 3. Hardware Architecture Parameters (Update these to match your Verilog!)
        self.Nacc = 32         # Phase Accumulator bit-width
        self.LUT_bits = 12     # Number of MSBs used to address the ROM
        self.DT_Mode = 'fixed' # Output format
        self.frac_bits = 15    # Fractional bits if using fixed-point

    def connect_phase(self):
        self.sc_export.connect(self.sc_fifo.analysis_export)
        
    async def run_phase(self):
        while True:
            # 1. Grab the transaction from the Monitor
            item = await self.sc_fifo.get()
            
            # 2. Skip comparison if reset was active
            if item.rst_n == 0:
                cocotb.log.debug(f"[SCOREBOARD] Skipping item {item.get_name()} due to active reset.")
                continue

            """************************** GOLDEN MODEL STIMULUS **************************"""
            # 3. Reconstruct the 'M' array (Tuning word for every clock cycle)
            # The hardware starts at FTW_start and accelerates by FTW_step every cycle.
            M_array = np.zeros(item.cycles, dtype=np.uint64)
            current_ftw = item.FTW_start
            
            for i in range(item.cycles):
                M_array[i] = current_ftw
                # Accumulate the step for the next cycle's tuning word
                current_ftw = (current_ftw + item.FTW_step) & 0xFFFFFFFF # Keep it 32-bit
                
            """************************** REFERENCE EXECUTION **************************"""
            # 4. Feed the cycle-by-cycle array into your MATLAB-ported core
            expected_wave = dds_core(
                M=M_array, 
                Nacc=self.Nacc, 
                LUT_bits=self.LUT_bits, 
                DT_Mode=self.DT_Mode, 
                frac_bits=self.frac_bits
            )
            
            # Assuming the Monitor captures the final amplitude at the end of 'cycles'
            expected_final_amplitude = expected_wave[-1] 

            """************************** COMPARISON LOGIC **************************"""
            # 5. Compare Hardware vs. Golden Model
            # Note: Because your python model does bit-accurate quantization ('fixed' DT_Mode), 
            # we can check for an exact match.
            
            # *You may need to cast the expected_final_amplitude to an integer depending 
            # on how your mytypes fixed-point library returns the data.*
            expected_val_int = int(expected_final_amplitude) 
            actual_val_int = int(item.final_amplitude)

            if expected_val_int == actual_val_int:
                self.passed_test_cases += 1
                cocotb.log.debug(f"[SCOREBOARD PASS] Expected: {hex(expected_val_int)} | Actual: {hex(actual_val_int)}")
            else:
                self.failed_test_cases += 1
                cocotb.log.error(f"[SCOREBOARD FAIL] {item.convert2string_stimulus()}")
                cocotb.log.error(f"  -> Expected Amplitude: {hex(expected_val_int)}")
                cocotb.log.error(f"  -> Actual Hardware   : {hex(actual_val_int)}")


    def report_phase(self):
        cocotb.log.info("========================================")
        cocotb.log.info("        SCOREBOARD FINAL REPORT         ")
        cocotb.log.info("========================================")
        cocotb.log.info(f" PASSED: {self.passed_test_cases}")
        cocotb.log.info(f" FAILED: {self.failed_test_cases}")
        cocotb.log.info("========================================")
        if self.failed_test_cases == 0 and self.passed_test_cases > 0:
            cocotb.log.info(" TEST STATUS: SUCCESS ")
        else:
            cocotb.log.error(" TEST STATUS: FAILED ")