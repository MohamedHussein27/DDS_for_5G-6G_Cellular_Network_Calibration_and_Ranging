from pyuvm import uvm_sequence
from fft_seq_item import fft_item
import math
from pyuvm import *


class SeqToneSweep(uvm_sequence):
    async def body(self):
        # Sweep through 5 different frequency bins
        for k in [10, 500, 1024, 2000, 4000]:
            for n in range(4096):
                item = fft_item("item")
                await self.start_item(item)
                angle = (2 * math.pi * k * n) / 4096
                item.rst_n = 1
                item.valid_in = 1
                # Scaling to 29490 (roughly 90% of full scale to avoid overflow)
                item.in_real = int(29490 * math.cos(angle))
                item.in_imag = int(29490 * math.sin(angle))
                await self.finish_item(item)

        # ========================================================
        # PIPELINE DRAIN (Executes ONCE after all 5 frames finish)
        # ========================================================
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0  # Drop valid!
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)