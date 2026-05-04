"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging
 
    Module: ifft_stress_test.py
 
    Description:
        Pushes the IFFT datapath to its numerical limits. Tests saturation
        handling at the positive and negative extremes of the Q8.8 word,
        rapid sign-reversal to stress carry/borrow propagation in the butterfly
        adder/subtractor, and a full-range ramp that exercises every twiddle
        ROM entry simultaneously with wide-spectrum stimulus.
 
    Sequences  (in order)
    ──────────────────────
        maxpositive_sequence   – every bin in_real = +32767, in_imag = 0
        maxnegative_sequence   – every bin in_real = -32768, in_imag = 0
        alternating_sequence   – +32767 / -32768 alternating every sample
        ramp_sequence          – in_real ramps 0 → 0.9 in Q8.8 over N samples
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
from ifft_maxpositive_seq import maxpositive_sequence
from ifft_maxnegative_seq import maxnegative_sequence
from ifft_alternating_seq import alternating_sequence
from ifft_sequences import ramp_sequence
from ifft_ramp_seq import ramp_sequence
 
 
@pyuvm.test()
class ifft_stress_test(ifft_base_test):
 
    def build_phase(self):
        super().build_phase()
        self.seq_rst     = reset_sequence.create("seq_rst")
        self.seq_maxpos  = maxpositive_sequence.create("seq_maxpos")
        self.seq_maxneg  = maxnegative_sequence.create("seq_maxneg")
        self.seq_alt     = alternating_sequence.create("seq_alt")
        self.seq_ramp    = ramp_sequence.create("seq_ramp")

        # test multi frames so we drive inputs and then go to the next sequence without waiting for the previous one to output the results. This is to test the IFFT block's ability to handle back-to-back sequences without idle cycles in between.
        self.seq_rst.multi_frame     = True
        self.seq_maxpos.multi_frame  = True
        self.seq_maxneg.multi_frame  = True
        self.seq_alt.multi_frame     = True
        self.seq_ramp.multi_frame    = True
        
 
    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")
 
        await self.generate_clock()
 
        self.logger.info("Running: reset_sequence")
        await self._run(self.seq_rst)
 
        self.logger.info("Running: maxpositive_sequence  (in_real = +32767 every sample)")
        await self._run(self.seq_maxpos)
 
        self.logger.info("Running: maxnegative_sequence  (in_real = -32768 every sample)")
        await self._run(self.seq_maxneg)
 
        self.logger.info("Running: alternating_sequence  (+32767 / -32768 alternating)")
        await self._run(self.seq_alt)
 
        self.logger.info("Running: ramp_sequence  (0 → 0.9 over 4096 samples)")
        await self._run(self.seq_ramp)
 
        self.drop_objection()
 
    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Positive and negative saturation")
        self.logger.info("                 Carry/borrow propagation under rapid sign reversal")
        self.logger.info("                 Full dynamic range and wide-spectrum twiddle ROM coverage")
        self.logger.info("---------------------------------------------------------")