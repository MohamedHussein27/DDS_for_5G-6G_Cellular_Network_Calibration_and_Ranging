import numpy as np
import random
from pyuvm import *
from ifft_item import ifft_item
from ifft_seq_utils import N, _frame, MAX_VAL, MIN_VAL
# ─────────────────────────────────────────────────────────────────────────────
# random_sequence  – random input values
# ─────────────────────────────────────────────────────────────────────────────
class random_sequence(uvm_sequence):
    def __init__(self, name="random_sequence", multi_frame=False):
        super().__init__(name)
        self.multi_frame = multi_frame
        self.latency = 12

    async def body(self):
        for _ in range(N):
            val = random.randint(MIN_VAL, MAX_VAL)
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in:
                    rst_n == 1 and valid_in == 1
            )
            self.seq_item.in_real = val
            self.seq_item.in_imag = 0
            await self.finish_item(self.seq_item)

        if not self.multi_frame:
            for _ in range(N + self.latency):
                self.seq_item = ifft_item.create("seq_item")
                await self.start_item(self.seq_item)
                self.seq_item.randomize_with(
                    lambda rst_n, valid_in: rst_n == 1 and valid_in == 0
                )
                self.seq_item.in_real = 0
                self.seq_item.in_imag = 0
                await self.finish_item(self.seq_item)