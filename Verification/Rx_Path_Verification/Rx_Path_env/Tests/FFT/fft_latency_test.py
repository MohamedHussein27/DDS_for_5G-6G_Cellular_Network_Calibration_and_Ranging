
import cocotb
from pyuvm import *
import pyuvm
from fft_base_test import base_test
from fft_reset_seq import FftResetSeq
from fft_latency_seq import SeqLatency



import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *

import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *


@pyuvm.test()
class test_fft_latency(base_test):
    """
    Test Case: Latency & Throughput Test
    Description: Drives exactly one full frame and measures pipeline depth.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def measure_pipeline_latency(self):
        """Bulletproof background task to act as a latency stopwatch."""
        dut = cocotb.top  # Safely grab the top module
        
        # 1. Wait for valid_in to go HIGH (using string comparison to avoid 'x' crashes)
        while (dut.valid_in.value) != 1:
            await RisingEdge(dut.clk)
        
        self.logger.info(" Latency Stopwatch Started!")
        latency_cycles = 0
        
        # 2. Count clock cycles until valid_out asserts
        while (dut.valid_out.value) != 1:
            await RisingEdge(dut.clk)
            latency_cycles += 1
            
        # 3. Report and Check
        self.logger.info(f" MEASURED PIPELINE LATENCY: {latency_cycles} clock cycles")
        
        if latency_cycles == 4095:
            self.logger.info(" Latency Check: PASS  (Exactly 4095 cycles)")
        else:
            self.logger.error(f" Latency Check: FAIL  (Expected 4095, got {latency_cycles})")

    async def run_phase(self):
        self.raise_objection()
        
        await self.generate_clock()
        await self.run_initial_setup()
        
        # Start the background stopwatch
        cocotb.start_soon(self.measure_pipeline_latency())
        
        # Start the Sequence
        self.logger.info("Starting Latency Test (4096 random samples)...")
        seq = SeqLatency("latency_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # Wait for the pipeline to flush
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await RisingEdge(cocotb.top.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()