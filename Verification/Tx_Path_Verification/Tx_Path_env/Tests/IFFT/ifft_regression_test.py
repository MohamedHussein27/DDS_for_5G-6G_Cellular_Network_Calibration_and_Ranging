"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: ifft_regression_test.py
 
    Description:
        Full regression test. Runs every sequence in a single simulation
        to achieve maximum functional and structural coverage in one pass.
        Use this in nightly regression. All other individual tests are
        subsets of this run.
 
    Sequences  (via fullregression_sequence — all 16 in order)
    ───────────────────────────────────────────────────────────
        reset → allzero → impulse → dc → singletone k=1 → singletone k=100
        → singletone k=2047 → twotone → maxpositive → maxnegative
        → alternating → nyquist → ramp → random → burst_512
        → resetmidframe → backtoback
"""
 
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging
 
from ifft_base_test  import ifft_base_test
from ifft_reset_seq  import reset_sequence
from ifft_regression_seq  import fullregression_sequence
 
 
@pyuvm.test()
class ifft_regression_test(ifft_base_test):
 
    def build_phase(self):
        super().build_phase()
        self.seq_rst        = reset_sequence.create("seq_rst")
        self.reg_seq        = fullregression_sequence.create("seq")

        # test multi frames so we drive inputs and then go to the next sequence without waiting for the previous one to output the results. This is to test the IFFT block's ability to handle back-to-back sequences without idle cycles in between.
        self.seq_rst.multi_frame    = False
        self.reg_seq.multi_frame        = False


 
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
 
        await self.generate_clock()
        self.logger.info("Running: reset_sequence")
        await self._run(self.seq_rst)

        self.logger.info("Running: fullregression_sequence  (all 16 sequences in order)")
        await self._run(self.reg_seq)
 
        self.drop_objection()
 
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Full regression — all corner cases in one run")
        self.logger.info("                 Maximum functional and structural coverage")
        self.logger.info("---------------------------------------------------------")