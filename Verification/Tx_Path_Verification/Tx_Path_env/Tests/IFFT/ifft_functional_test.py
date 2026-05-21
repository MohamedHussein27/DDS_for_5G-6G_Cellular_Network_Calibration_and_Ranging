"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: ifft_functional_test.py
 
    Description:
        Verifies the core IFFT mathematical correctness across the most
        fundamental input patterns. Each sequence produces a well-defined,
        analytically predictable output that the scoreboard checks against
        the numpy golden model.
 
    Sequences  (in order)
    ──────────────────────
        allzero_sequence        – all bins 0+0j          → all-zero output
        impulse_sequence        – bin 0 = 1+0j only      → constant 1/N for all n
        dc_sequence             – all bins = 1+0j        → spike at n=0, zeros elsewhere
        singletone_sequence k=1 – unit phasor at bin 1   → pure sinusoid at f=1/N
        twotone_sequence        – bins k=10 and k=200    → two-sinusoid superposition
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
from ifft_allzero_seq import allzero_sequence
from ifft_impulse_seq import impulse_sequence
from ifft_dc_seq    import dc_sequence
from ifft_singletone_seq import singletone_sequence
from ifft_twotone_seq import twotone_sequence
 
@pyuvm.test()
class ifft_functional_test(ifft_base_test):
 
    def build_phase(self):
        super().build_phase()
        # Create the objects using the strict pyuvm factory
        self.seq_rst       = reset_sequence.create("seq_rst")
        self.seq_allzero   = allzero_sequence.create("seq_allzero")
        self.seq_impulse   = impulse_sequence.create("seq_impulse")
        self.seq_dc        = dc_sequence.create("seq_dc")
        self.seq_tone_k1   = singletone_sequence.create("seq_tone_k1")
        self.seq_twotone   = twotone_sequence.create("seq_twotone")

        # test multi frames so we drive inputs and then go to the next sequence without waiting for the previous one to output the results. This is to test the IFFT block's ability to handle back-to-back sequences without idle cycles in between.
        self.seq_rst.multi_frame     = False
        self.seq_allzero.multi_frame = False
        self.seq_impulse.multi_frame = False
        self.seq_dc.multi_frame      = False
        self.seq_tone_k1.multi_frame = False
        self.seq_twotone.multi_frame = False
        
        # Other existing assignments
        self.seq_tone_k1.k = 1
 
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
 
        await self.generate_clock()

        self.logger.info("Running: reset_sequence") # always run reset first to ensure a known starting state for the IFFT block
        await self._run(self.seq_rst)

        self.logger.info("Running: allzero_sequence  (zero input → zero output)")
        await self._run(self.seq_allzero)

        self.logger.info("Running: impulse_sequence  (bin 0 → constant 1/N output)")
        await self._run(self.seq_impulse)

        self.logger.info("Running: dc_sequence  (all bins = 1 → spike at n=0)")
        await self._run(self.seq_dc)
        
        self.logger.info("Running: singletone_sequence  k=1  (pure sinusoid at f=1/N)")
        await self._run(self.seq_tone_k1)
 
        self.logger.info("Running: twotone_sequence  (bins k=10, k=200 → two sinusoids)")
        await self._run(self.seq_twotone)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Core IFFT math — zero, impulse, DC, single-tone,")
        self.logger.info("                 two-tone linearity")
        self.logger.info("---------------------------------------------------------")