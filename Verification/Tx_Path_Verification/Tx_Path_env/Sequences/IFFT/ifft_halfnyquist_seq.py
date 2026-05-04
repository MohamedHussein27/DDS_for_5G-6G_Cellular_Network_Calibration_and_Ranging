import numpy as np
from pyuvm import *
from ifft_item import ifft_item
from ifft_seq_utils import N, _frame, MAX_VAL, MIN_VAL
# ─────────────────────────────────────────────────────────────────────────────
# halfnyquist_sequence  – bin N/2 = 2048 → alternating ±1/N output
# ─────────────────────────────────────────────────────────────────────────────
class halfnyquist_sequence(uvm_sequence):
    def __init__(self, name="halfnyquist_sequence", multi_frame=False):
        super().__init__(name)
        self.multi_frame = multi_frame

    async def body(self):
        X = np.zeros(N, dtype=complex)
        X[N // 2] = 1.0 + 0j
        pairs = _frame(X)

        for re, im in pairs:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in:
                    rst_n == 1 and valid_in == 1 
            )
            self.seq_item.in_real = re
            self.seq_item.in_imag = im  
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