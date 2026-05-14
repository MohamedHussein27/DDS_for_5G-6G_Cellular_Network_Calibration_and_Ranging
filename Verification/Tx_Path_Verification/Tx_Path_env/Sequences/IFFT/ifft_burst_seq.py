import numpy as np
from pyuvm import *
from ifft_item import ifft_item
from ifft_seq_utils import N, _frame, MAX_VAL, MIN_VAL
# ─────────────────────────────────────────────────────────────────────────────
# burst_sequence  – burst of activity then idle
# ─────────────────────────────────────────────────────────────────────────────
class burst_sequence(uvm_sequence):
    def __init__(self, name="burst_sequence", burst_length=10, multi_frame=False):
        super().__init__(name)
        self.burst_length = burst_length
        self.multi_frame = multi_frame
        self.latency = 12

    async def body(self):
        # Burst of valid data
        for i in range(min(self.burst_length, N)):
            val = MAX_VAL if i % 2 == 0 else MIN_VAL
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in:
                    rst_n == 1 and valid_in == 1
            )
            self.seq_item.in_real = val
            self.seq_item.in_imag = val
            await self.finish_item(self.seq_item)

        # Remaining idle cycles
        for _ in range(N - self.burst_length):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            self.seq_item.randomize_with(
                lambda rst_n, valid_in: rst_n == 1 and valid_in == 0
            )
            self.seq_item.in_real = 0
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