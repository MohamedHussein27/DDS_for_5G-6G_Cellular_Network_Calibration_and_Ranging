import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
import logging
from ifft_item import ifft_item
from ifft_seq_utils import *
from ifft_reset_seq import reset_sequence
from ifft_allzero_seq import allzero_sequence
from ifft_impulse_seq import impulse_sequence
from ifft_dc_seq import dc_sequence
from ifft_singletone_seq import singletone_sequence
from ifft_twotone_seq import twotone_sequence
from ifft_maxpositive_seq import maxpositive_sequence
from ifft_maxnegative_seq import maxnegative_sequence   
from ifft_alternating_seq import alternating_sequence
from ifft_halfnyquist_seq import halfnyquist_sequence
from ifft_ramp_seq import ramp_sequence
from ifft_random_seq import random_sequence
from ifft_burst_seq import burst_sequence
from ifft_resetmidframe_seq import resetmidframe_sequence
from ifft_backtoback_seq import backtoback_sequence


# ─────────────────────────────────────────────────────────────────────────────
# 16. fullregression_sequence  – runs every sequence in order
# ─────────────────────────────────────────────────────────────────────────────
class fullregression_sequence(uvm_sequence):
    def __init__(self, name="fullregression_sequence", multi_frame=True):
        super().__init__(name)
        self.multi_frame = multi_frame


    async def body(self):
        sequences = [
            reset_sequence        ("reset", multi_frame=self.multi_frame),
            allzero_sequence      ("allzero", multi_frame=self.multi_frame),
            impulse_sequence      ("impulse", multi_frame=self.multi_frame),
            dc_sequence           ("dc", multi_frame=self.multi_frame),
            singletone_sequence   ("tone_k1",    k=1, amplitude=0.625 , multi_frame=self.multi_frame),
            singletone_sequence   ("tone_k100",  k=100, amplitude=0.125, multi_frame=self.multi_frame),
            singletone_sequence   ("tone_k2047", k=2047, amplitude=0.0625, multi_frame=self.multi_frame),
            twotone_sequence      ("twotone", multi_frame=self.multi_frame),
            maxpositive_sequence  ("maxpos", multi_frame=self.multi_frame),
            maxnegative_sequence  ("maxneg", multi_frame=self.multi_frame),
            alternating_sequence  ("alternating", multi_frame=self.multi_frame),
            halfnyquist_sequence  ("nyquist", multi_frame=self.multi_frame),
            ramp_sequence         ("ramp", multi_frame=self.multi_frame),
            random_sequence       ("random", multi_frame=self.multi_frame),
            burst_sequence        ("burst_512",  burst_length=512, multi_frame=self.multi_frame),
            #backtoback_sequence   ("b2b"),
        ]

        for seq in sequences:
            print(f"[fullregression_sequence] ── Starting {seq.get_type_name()} ──")
            await seq.start(self.sequencer)