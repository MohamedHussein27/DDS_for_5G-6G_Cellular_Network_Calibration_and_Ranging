"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_singletone_test.py
    Description: Verifies a single active bin produces a pure sinusoid.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
from ifft_base_test  import ifft_base_test
from ifft_sequences  import reset_sequence, singletone_sequence

@pyuvm.test()
class ifft_singletone_test(ifft_base_test):

    def build_phase(self):
        super().build_phase()
        self.seq_rst     = reset_sequence.create("seq_rst")
        self.seq_tone_k1 = singletone_sequence.create("seq_tone_k1")

        self.seq_rst.multi_frame     = False
        self.seq_tone_k1.multi_frame = False
        self.seq_tone_k1.k           = 1

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        await self.generate_clock()

        self.logger.info("Running: reset_sequence")
        await self._run(self.seq_rst)

        self.logger.info("Running: singletone_sequence  k=1  (pure sinusoid at f=1/N)")
        await self._run(self.seq_tone_k1)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Core IFFT math — Single Tone")
        self.logger.info("---------------------------------------------------------")