import numpy as np
from pyuvm import *
from ifft_item import ifft_item
from ifft_seq_utils import N, _frame

class twotone_sequence(uvm_sequence):
    def __init__(self, name="twotone_sequence", k1=10, k2=200, a1=0.5, a2=0.5, multi_frame=False):
        super().__init__(name)
        self.k1, self.k2 = k1, k2
        self.a1, self.a2 = a1, a2
        self.multi_frame = multi_frame

    async def body(self):
        X = np.zeros(N, dtype=complex)
        X[self.k1] = self.a1
        X[self.k2] = self.a2
        pairs = _frame(X)

        # ── First N cycles: drive the input frame ──
        for re, im in pairs:
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            
            self.seq_item.randomize_with(
                lambda rst_n, valid_in: rst_n == 1 and valid_in == 1
            )
            self.seq_item.in_real = re
            self.seq_item.in_imag = im
            
            await self.finish_item(self.seq_item)

        # ── Second N cycles: idle while output comes out (Drain time) ──
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