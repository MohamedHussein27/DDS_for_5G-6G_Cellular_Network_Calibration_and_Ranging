"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_twotone_test.py
    Description: Verifies two active bins produce the correct superposition.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
from ifft_base_test  import ifft_base_test
from ifft_sequences  import reset_sequence, twotone_sequence

@pyuvm.test()
class ifft_twotone_test(ifft_base_test):

    def build_phase(self):
        super().build_phase()
        self.seq_rst     = reset_sequence.create("seq_rst")
        self.seq_twotone = twotone_sequence.create("seq_twotone")

        self.seq_rst.multi_frame     = False
        self.seq_twotone.multi_frame = False

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        await self.generate_clock()

        self.logger.info("Running: reset_sequence")
        await self._run(self.seq_rst)

        self.logger.info("Running: twotone_sequence  (bins k=10, k=200 → two sinusoids)")
        await self._run(self.seq_twotone)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Core IFFT math — Two Tone Linearity")
        self.logger.info("---------------------------------------------------------")