import pyuvm
from pyuvm import *
from dds_seq_item import *
import cocotb
import numpy as np

# 1. Import the Golden Model
from dds_golden_model import DDSGoldenModel 

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
        # 14-bit addressing (16384) mapped to 16-bit phase (65536) for a quarter wave
        sine_float = np.sin(2.0 * np.pi * addresses / 65536.0)
        self.rom_data = np.round(sine_float * 127).astype(np.int32)

        # 2. Instantiate the Golden Model
        self.golden_model = DDSGoldenModel(rom_data=self.rom_data)

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
                    # 3. Call the model once per chirp sequence
                    self.expected_wave = self.golden_model.generate_chirp(
                        ftw_start=item.FTW_start, 
                        ftw_step=item.FTW_step, 
                        cycles=item.cycles
                    )

                """********* COMPARISON LOGIC *********"""
                # No more ignoring cycles! valid_out guarantees this is perfect data.
                math_idx = item.sample_index - 1
                
                if math_idx < len(self.expected_wave):
                    expected_val_int = int(self.expected_wave[math_idx])
                    actual_val_int = int(item.final_amplitude)

                    if expected_val_int == actual_val_int:
                        self.passed_test_cases += 1
                        cocotb.log.info(f"[SCOREBOARD PASS] Stream Cycle {item.sample_index} Match! Expected: {hex(expected_val_int & 0xFF)} | Actual: {hex(actual_val_int & 0xFF)}")
                    else:
                        self.failed_test_cases += 1
                        cocotb.log.error(f"[SCOREBOARD FAIL] Mismatch at Stream Cycle {item.sample_index}!")
                        cocotb.log.error(f"  -> {item.convert2string_stimulus()}")
                        cocotb.log.error(f"  -> Expected Amplitude: {hex(expected_val_int & 0xFF)}")
                        cocotb.log.error(f"  -> Actual Hardware   : {hex(actual_val_int & 0xFF)}")

    def report_phase(self):
        cocotb.log.info("========================================")
        cocotb.log.info("       STREAMING SCOREBOARD REPORT      ")
        cocotb.log.info("========================================")
        cocotb.log.info(f" SAMPLES PASSED: {self.passed_test_cases}")
        cocotb.log.info(f" SAMPLES FAILED: {self.failed_test_cases}")
        cocotb.log.info("========================================")
        if self.failed_test_cases == 0 and self.passed_test_cases > 0:
            cocotb.log.info(" TEST STATUS: SUCCESS ")
        else:
            cocotb.log.error(" TEST STATUS: FAILED ")