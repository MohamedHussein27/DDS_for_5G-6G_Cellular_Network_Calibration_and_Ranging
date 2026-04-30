import cocotb
from pyuvm import *
import pyuvm
from base_test import base_test
from fft_reset_seq import FftResetSeq
#from fft_zero_seq import SeqZero
from fft_impulse_seq import SeqImpulse
from fft_ten_sample_seq import fft_ten_sample_seq
from fft_multiframe_seq import fft_two_frame_seq
from fft_dc_seq import SeqDc
from fft_latency_seq import SeqLatency
from fft_neg_dc_seq import SeqNegDc
from fft_fullscale_seq import SeqFullScale
# --- Sequence Imports ---
#import fft_random_seq
#from fft_latency_seq import SeqLatency  # <-- Assuming your SeqLatency class is in this file


# ==============================================================================
# TEST 2: Latency & Throughput Test
# Run via: make -f runs.mk run_latency
# =============================================================================

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
        # TEST 3: Mid-stream Reset Recovery Test
# ==============================================================================
@pyuvm.test()
class test_fft_recovery(base_test):
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        await self.generate_clock()
        await self.run_initial_setup() # This runs the safe boot-up reset
        
        self.logger.info("Starting Reset Recovery Stress Test...")
        seq = FftResetSeq("recovery_seq") # This runs YOUR sequence!
        await seq.start(self.env.fft_agt.sqr)
        
        for _ in range(4500):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.drop_objection()


        # ==============================================================================
# TEST 4: Zero Signal (Silence) Test
# Run via: make -f Makefile.runs run_zero
# ==============================================================================
@pyuvm.test()
class test_fft_zero(base_test):
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the Zero Sequence
        self.logger.info("Starting Zero-Input (Silence) Test...")
        seq = SeqZero("zero_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the pipeline to flush (N-1 cycles)
        # We need to see if the FFT actually outputs 0s after the latency
        for _ in range(4500):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.drop_objection()
        
# TEST 5: Impulse Response (Kronecker Delta) Test
# Run via: make -f Makefile.runs run_impulse
# ==============================================================================
@pyuvm.test()
class test_fft_impulse(base_test):
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # Start Hardware Clock & Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # Start the Impulse Sequence
        self.logger.info("Starting Impulse Response (Delta) Test...")
        seq = SeqImpulse("impulse_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # Wait for the pipeline to flush
        # This is critical to see the frequency response results in the scoreboard
        for _ in range(4500):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.drop_objection()


@pyuvm.test()
class test_fft_ten_samples(base_test):
    """
    Test Case: Ten Sample Input with Zero-Padding
    Description: Drives 10 custom real-valued samples followed by zeros 
    to verify the FFT's transient response and frequency spreading.
    """
    def build_phase(self):
        # Set the mode so the Environment knows we are verifying the FFT block
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Global Reset from base_test
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the Ten-Sample Sequence
        self.logger.info("Starting Ten-Sample Input (Zero-Padded) Test...")
        
        # Ensure the sequence name matches your sequence file class name
        seq = fft_ten_sample_seq("ten_sample_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the pipeline to flush
        # For a 4096-point FFT, we wait >4096 cycles to ensure the 
        # result for our 10 samples propagates through all internal stages.
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(4500):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()

@pyuvm.test()
class test_fft_multiframe(base_test):
    """
    Test Case: Multi-Frame (Two Consecutive Frames)
    Description: Drives two back-to-back 4096-point frames to verify 
    the pipeline can handle continuous data without dropping samples 
    between frames.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Global Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the Multi-Frame Sequence
        self.logger.info("Starting Multi-Frame (Continuous) Test...")
        
        # Ensure this matches the class name you used in your sequence file!
        seq = fft_two_frame_seq("two_frame_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the pipeline to drain
        # TWO frames = 8192 valid_in cycles. 
        # We must wait >8192 cycles + the RTL pipeline latency.
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(10000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()


@pyuvm.test()
class test_fft_dc(base_test):
    """
    Test Case: DC Input Test
    Description: Drives a constant signal (29490) across all 4096 samples 
    to verify the DC bin (Bin 0) accumulates correctly and no other 
    frequency bins show phantom energy.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Global Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the DC Sequence
        self.logger.info("Starting DC Input Test (Constant 29490)...")
        
        # Ensure this matches the class name in your sequence file
        seq = SeqDc("dc_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the pipeline to flush
        # 4096 cycles for data to enter + RTL pipeline latency
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()


@pyuvm.test()
class test_fft_neg_dc(base_test):
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        await self.generate_clock()
        await self.run_initial_setup()
        
        self.logger.info("Starting Negative DC Test (-32768)...")
        seq = SeqNegDc("neg_dc_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()


@pyuvm.test()
class test_fft_full_scale(base_test):
    """
    Test Case: Full-Scale Nyquist Stress Test
    Description: Pushes alternating max/min values (+32767, -32768) 
    to test the highest possible frequency (Nyquist) and verify 
    the internal multipliers and adders do not overflow.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        await self.generate_clock()
        await self.run_initial_setup()
        
        self.logger.info("Starting Full-Scale Nyquist Test...")
        
        # Make sure this matches the class name in your sequence file!
        seq = SeqFullScale("full_scale_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        self.logger.info("Sequence finished. Waiting for pipeline flush...")
        for _ in range(5000):
            await cocotb.triggers.RisingEdge(self.dut.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()

