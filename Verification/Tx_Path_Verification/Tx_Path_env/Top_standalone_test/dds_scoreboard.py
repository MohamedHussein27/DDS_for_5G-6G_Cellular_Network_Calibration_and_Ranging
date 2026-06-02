
import pyuvm
from pyuvm import *
from dds_seq_item import *
import cocotb
import numpy as np

class dds_scoreboard(uvm_scoreboard):
    def build_phase(self):
        super().build_phase()
        self.passed_test_cases = 0
        self.failed_test_cases = 0
        
        self.sb_export = uvm_analysis_export("sb_export", self)
        self.sb_fifo = uvm_tlm_analysis_fifo("sb_fifo", self)
        self.sb_export = self.sb_fifo.analysis_export
        
        self.Nacc = 32         
        self.LUT_bits = 12     

        cocotb.log.info("[SCOREBOARD] Generating mathematical ROM natively...")
        addresses = np.arange(16384)
        sine_float = np.sin(2.0 * np.pi * addresses / 65536.0)
        self.rom_data = np.round(sine_float * 127).astype(np.int32)

    async def run_phase(self):
        while True:
            item = await self.sb_fifo.get()
            
            # 1. Verify Reset State Instantly
            if item.rst_n == 0:
                if item.final_amplitude == 0:
                    self.passed_test_cases += 1
                else:
                    self.failed_test_cases += 1
                    cocotb.log.error(f"[SCOREBOARD FAIL] Hardware Reset Failed! Expected 0x0, got {hex(item.final_amplitude & 0xFF)}")
                continue

            if item.cycles <= 1:
                continue

            # 2. Verify Streaming Data
            if item.valid_out == 1:
                """********* NATIVE BIT-TRUE PYTHON MODEL *********"""
                if item.sample_index == 1:
                    # Start 'n' at 1, Verilog skips 0 and jumps straight to FTW_start
                    n = np.arange(1, item.cycles + 1, dtype=np.uint64)
                    FTW_start = np.uint64(item.FTW_start)
                    FTW_step = np.uint64(item.FTW_step)
                    
                    accumulation_term = (n * (n - np.uint64(1))) // np.uint64(2)
                    ideal_discrete_phase_word = (FTW_start * n + FTW_step * accumulation_term) & 0xFFFFFFFF
                    
                    truncated_phase_16b = ideal_discrete_phase_word >> 16
                    quadrant = truncated_phase_16b >> 14
                    addr = truncated_phase_16b & 0x3FFF
                    
                    mapped_addr = np.where((quadrant == 1) | (quadrant == 3), 16383 - addr, addr)
                    neg_flag = np.where((quadrant == 2) | (quadrant == 3), 1, 0)
                    
                    lut_amplitude = self.rom_data[mapped_addr]
                    self.expected_wave = np.where(neg_flag == 1, -lut_amplitude, lut_amplitude)

                """********* COMPARISON LOGIC *********"""
                # No more ignoring cycles! valid_out guarantees this is perfect data.
                math_idx = item.sample_index - 1
                
                if math_idx < len(self.expected_wave):
                    expected_val_int = int(self.expected_wave[math_idx])
                    actual_val_int = int(item.final_amplitude)

                    if expected_val_int == actual_val_int:
                        self.passed_test_cases += 1
                        cocotb.log.debug(f"[SCOREBOARD PASS] Stream Cycle {item.sample_index} Match! Expected: {hex(expected_val_int & 0xFF)} | Actual: {hex(actual_val_int & 0xFF)}")
                    else:
                        self.failed_test_cases += 1
                        cocotb.log.error(f"[SCOREBOARD FAIL] Mismatch at Stream Cycle {item.sample_index}!")
                        cocotb.log.error(f"  -> {item.convert2string_stimulus()}")
                        cocotb.log.error(f"  -> Expected Amplitude: {hex(expected_val_int & 0xFF)}")
                        cocotb.log.error(f"  -> Actual Hardware   : {hex(actual_val_int & 0xFF)}")

    def report_phase(self):
        cocotb.log.info("========================================")
        cocotb.log.info("      DDS SCOREBOARD REPORT       ")
        cocotb.log.info("========================================")
        cocotb.log.info(f" SAMPLES PASSED: {self.passed_test_cases}")
        cocotb.log.info(f" SAMPLES FAILED: {self.failed_test_cases}")
        cocotb.log.info("========================================")
        if self.failed_test_cases == 0 and self.passed_test_cases > 0:
            cocotb.log.info(" TEST STATUS: SUCCESS ")
        else:
            cocotb.log.error(" TEST STATUS: FAILED ")