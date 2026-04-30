import cocotb
from cocotb.triggers import RisingEdge
import pyuvm
from pyuvm import *

# Pull in the base test and your new Tone Sweep sequence
from fft_base_test import base_test
from fft_tone_sweep_seq import SeqToneSweep
from fft_reset_seq import FftResetSeq



@pyuvm.test()
class test_fft_tone_sweep(base_test):
    """
    Test Case: Single-Tone Frequency Sweep
    Description: Drives 5 consecutive frames of perfect sine waves at 
    different frequency bins to verify spectral isolation and peak accuracy.
    """
    def build_phase(self):
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Start Hardware Clock & Global Reset
        await self.generate_clock()
        await self.run_initial_setup()
        
        # 2. Start the Sequence
        self.logger.info("Starting Tone Sweep Test (5 Frames)...")
        seq = SeqToneSweep("tone_sweep_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Wait for the hardware pipeline to process the final frame
        self.logger.info("Sequence finished. Waiting for final pipeline flush...")
        for _ in range(5000):
            await RisingEdge(cocotb.top.clk)
            
        self.logger.info("Test Complete. Dropping objection.")
        self.drop_objection()