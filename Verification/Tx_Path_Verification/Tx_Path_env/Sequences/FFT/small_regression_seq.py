from pyuvm import uvm_sequence
import cocotb
from pyuvm import *
# Sequence Imports
from fft_impulse_seq import SeqImpulse
from fft_dc_seq import SeqDc
from fft_ten_sample_seq import fft_ten_sample_seq


class SeqSmallRegression(uvm_sequence):
    """
    Executes a mini-regression of 3 specific FFT sequences back-to-back.
    """

    async def body(self):

        cocotb.log.info("=== STARTING SMALL REGRESSION ===")

        sequence_list = [
            (SeqImpulse, "Impulse Test"),
            (SeqDc, "DC Test"),
            (fft_ten_sample_seq, "Ten Sample Test")
        ]

        for SeqClass, test_name in sequence_list:
            cocotb.log.info("======================================")
            cocotb.log.info(f" RUNNING: {test_name}")
            cocotb.log.info("======================================")

            seq = SeqClass(f"seq_{test_name.replace(' ', '_')}")
            await seq.start(self.sequencer)   #  FIXED (no parentheses)

        cocotb.log.info("=== SMALL REGRESSION COMPLETE ===")

   
