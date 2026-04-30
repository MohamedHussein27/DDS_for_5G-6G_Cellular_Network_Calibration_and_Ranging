from pyuvm import *
from fft_seq_item import fft_item

class SeqDc(uvm_sequence):
    async def body(self):
        for _ in range(4096):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = 29490
            item.in_imag = 0
            await self.finish_item(item)