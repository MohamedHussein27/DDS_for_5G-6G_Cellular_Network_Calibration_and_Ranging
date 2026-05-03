"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_allzero_test.py
    Description: Verifies that an all-zero input frame produces an all-zero output.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
from ifft_base_test  import ifft_base_test
from ifft_sequences  import reset_sequence, allzero_sequence

@pyuvm.test()
class ifft_allzero_test(ifft_base_test):

    def build_phase(self):
        super().build_phase()
        self.seq_rst     = reset_sequence.create("seq_rst")
        self.seq_allzero = allzero_sequence.create("seq_allzero")

        self.seq_rst.multi_frame     = False
        self.seq_allzero.multi_frame = False

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        await self.generate_clock()

        self.logger.info("Running: reset_sequence")
        await self._run(self.seq_rst)

        self.logger.info("Running: allzero_sequence  (zero input → zero output)")
        await self._run(self.seq_allzero)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Core IFFT math — All Zeroes")
        self.logger.info("---------------------------------------------------------")