"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: ifft_frequency_test.py
 
    Description:
        Sweeps energy across the full frequency axis of the IFFT to verify
        twiddle-factor ROM accuracy, stage routing, and bin-address decoding
        from DC (k=0) to Nyquist (k=N/2). Includes three single-tone bins
        at low, mid, and high frequencies, a two-tone linearity check, and
        the Nyquist corner case where the twiddle factor W = -1.
 
    Sequences  (in order)
    ──────────────────────
        singletone_sequence k=1    – low frequency  (tests lower twiddle ROM entries)
        singletone_sequence k=100  – mid frequency
        singletone_sequence k=2047 – high frequency (just below Nyquist)
        twotone_sequence           – k=10 and k=200 simultaneously (linearity)
        halfnyquist_sequence       – k=2048 Nyquist bin (W = -1 corner case)
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
from ifft_singletone_seq import singletone_sequence
from ifft_twotone_seq import twotone_sequence
from ifft_halfnyquist_seq import halfnyquist_sequence

 
 
@pyuvm.test()
class ifft_frequency_test(ifft_base_test):
 
    def build_phase(self):
        super().build_phase()
        self.seq_rst      = reset_sequence.create("seq_rst")
        self.seq_k1      = singletone_sequence.create("seq_k1")
        self.seq_k1.k    = 1
 
        self.seq_k100    = singletone_sequence.create("seq_k100")
        self.seq_k100.k  = 100
 
        self.seq_k2047   = singletone_sequence.create("seq_k2047")
        self.seq_k2047.k = 2047
 
        self.seq_twotone = twotone_sequence.create("seq_twotone")
        self.seq_nyquist = halfnyquist_sequence.create("seq_nyquist")

        # test multi frames so we drive inputs and then go to the next sequence without waiting for the previous one to output the results. This is to test the IFFT block's ability to handle back-to-back sequences without idle cycles in between.
        self.seq_rst.multi_frame        = True
        self.seq_k100.multi_frame       = True
        self.seq_k2047.multi_frame      = True
        self.seq_twotone.multi_frame    = True
        self.seq_nyquist.multi_frame    = True
        
        
 
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
 
        await self.generate_clock()

        self.logger.info("Running: reset_sequence")
        await self._run(self.seq_rst)
 
        self.logger.info("Running: singletone_sequence  k=1  (low frequency)")
        await self._run(self.seq_k1)
 
        self.logger.info("Running: singletone_sequence  k=100  (mid frequency)")
        await self._run(self.seq_k100)
 
        self.logger.info("Running: singletone_sequence  k=2047  (high frequency, just below Nyquist)")
        await self._run(self.seq_k2047)
 
        self.logger.info("Running: twotone_sequence  (k=10 and k=200, linearity check)")
        await self._run(self.seq_twotone)
 
        self.logger.info("Running: halfnyquist_sequence  (k=2048, W=-1 corner case)")
        await self._run(self.seq_nyquist)
 
        self.drop_objection()
 
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Twiddle-ROM accuracy at bins 1, 100, 2047")
        self.logger.info("                 Two-tone linearity, Nyquist bin corner case")
        self.logger.info("---------------------------------------------------------")