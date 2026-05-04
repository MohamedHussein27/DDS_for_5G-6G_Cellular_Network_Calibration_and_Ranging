import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *
import pyuvm

# Import your base test and the new sequence
from fft_base_test import base_test
from fft_coverage_closure_seq import FftCoverageClosureSeq

@pyuvm.test()
class test_fft_coverage_closure(base_test):
    """
    Test Case: Coverage Closure Test
    Description: Executes the 25-permutation matrix and pure DC 
    saturation frames (Max Positive and Min Negative) to hit the 
    absolute extreme boundary conditions and close all functional 
    coverage holes for the 16-bit FFT datapath.
    """
    def build_phase(self):
        # Set the environment to run the standalone FFT
        ConfigDB().set(self, "*", "VERIF_MODE", "FFT")
        super().build_phase()

    async def run_phase(self):
        self.raise_objection()
        
        # 1. Startup Routine
        await self.generate_clock()
        await self.run_initial_setup()
        
        self.logger.info("Starting Coverage Closure Test...")
        
        # 2. Execute the targeted sequence
        seq = FftCoverageClosureSeq("coverage_closure_seq")
        await seq.start(self.env.fft_agt.sqr)
        
        # 3. Final Pipeline Drain
        # (The sequence already pushes 100 idle cycles, but we add 
        # a standard wait here to ensure the monitor captures everything)
        self.logger.info("Sequence finished. Waiting for final pipeline flush...")
        for _ in range(5000):
            await RisingEdge(self.dut.clk)
            
        self.logger.info("Coverage Closure Test Complete. Dropping objection.")
        self.drop_objection()