"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_maxnegative_test.py
    Description: Verifies maximum negative input value (-32768) behavior.
"""
import cocotb
from cocotb.triggers import *
from cocotb.clock    import Clock
from cocotb_coverage.crv import *
from pyuvm import *
import pyuvm
from ifft_base_test  import ifft_base_test
from ifft_reset_seq  import reset_sequence
from ifft_maxnegative_seq import maxnegative_sequence

@pyuvm.test()
class ifft_maxnegative_test(ifft_base_test):

    def build_phase(self):
        super().build_phase()
        self.seq_rst           = reset_sequence.create("seq_rst")
        self.seq_maxnegative   = maxnegative_sequence.create("seq_maxnegative")

        self.seq_rst.multi_frame         = False
        self.seq_maxnegative.multi_frame = False

    async def run_phase(self):
        self.raise_objection()
        self.logger.info(f"================ Start of {self.get_type_name()} ================")

        await self.generate_clock()

        self.logger.info("Running: reset_sequence")
        await self._run(self.seq_rst)

        self.logger.info("Running: maxnegative_sequence (all samples = -32768)")
        await self._run(self.seq_maxnegative)

        self.drop_objection()

    def report_phase(self):
        self.logger.info("---------------------------------------------------------")
        self.logger.info(f" [TEST REPORT] : {self.get_type_name()}")
        self.logger.info(" [TARGETS]     : Core IFFT math — Maximum Negative Input")
        self.logger.info("---------------------------------------------------------")