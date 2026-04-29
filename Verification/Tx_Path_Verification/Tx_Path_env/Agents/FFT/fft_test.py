import cocotb
from pyuvm import *
import pyuvm
from base_test import base_test
from fft_reset_seq import FftResetSeq
#from fft_zero_seq import SeqZero
from fft_impulse_seq import SeqImpulse
from fft_ten_sample_seq import fft_ten_sample_seq
# --- Sequence Imports ---
#import fft_random_seq
#from fft_latency_seq import SeqLatency  # <-- Assuming your SeqLatency class is in this file


# ==============================================================================
# TEST 2: Latency & Throughput Test
# Run via: make -f runs.mk run_latency
# ==============================================================================
@pyuvm.test()
class test_fft_latency(base_test):
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # Start Hardware Clock & Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # Start the Latency Sequence (4096 back-to-back valid inputs)
        self.logger.info("Starting Latency Sequence...")
        seq = SeqLatency("latency_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        
            
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