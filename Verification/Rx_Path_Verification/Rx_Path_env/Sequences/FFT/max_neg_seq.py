import cocotb
from pyuvm import uvm_sequence
from fft_seq_item import fft_item

from pyuvm import *


class FftUnscaledMinSeq(uvm_sequence):
    """
    User-designed sequence for an UNSCALED 4096-point FFT.
    Injects a continuous DC frame of -8 to perfectly sum to -32768 
    at the output DC bin, closing the out_real_min_neg coverage hole.
    """
    async def body(self):
        N = 4096
        
        # ==========================================
        # PHASE 1: The "-8" DC Frame
        # ==========================================
        cocotb.log.info("Injecting -8 DC Frame to hit -32768 output...")
        for i in range(N):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = -8       # <--- Your exact solution!
            item.in_imag = 0        # Keep imaginary quiet
            await self.finish_item(item)

        # ==========================================
        # PHASE 2: PIPELINE DRAIN
        # ==========================================
        cocotb.log.info("Flushing Pipeline...")
        for i in range(N + 100):
            idle_item = fft_item("idle")
            await self.start_item(idle_item)
            idle_item.rst_n = 1
            idle_item.valid_in = 0
            idle_item.in_real = 0
            idle_item.in_imag = 0
            await self.finish_item(idle_item)