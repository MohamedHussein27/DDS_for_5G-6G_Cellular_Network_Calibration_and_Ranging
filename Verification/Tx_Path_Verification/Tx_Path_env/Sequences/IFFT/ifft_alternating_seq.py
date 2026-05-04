import numpy as np
from pyuvm import *
from ifft_item import ifft_item
from ifft_seq_utils import N, _frame, MAX_VAL, MIN_VAL
# ─────────────────────────────────────────────────────────────────────────────
# alternating_sequence  – +MAX / -MAX every other sample
# ─────────────────────────────────────────────────────────────────────────────
class alternating_sequence(uvm_sequence):
    def __init__(self, name="alternating_sequence", multi_frame=False):
        super().__init__(name)
        self.multi_frame = multi_frame

    async def body(self):
        for i in range(N):
            val = MAX_VAL if i % 2 == 0 else MIN_VAL
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
            for _ in range(N):
                self.seq_item = ifft_item.create("seq_item")
                await self.start_item(self.seq_item)
                self.seq_item.randomize_with(
                    lambda rst_n, valid_in: rst_n == 1 and valid_in == 0
                )
                self.seq_item.in_real = 0
                self.seq_item.in_imag = 0
                await self.finish_item(self.seq_item)