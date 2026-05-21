import random
from pyuvm import *
from fft_seq_item import fft_item

class SeqLatency(uvm_sequence):
    async def body(self):
        # 1. Send the 4096 sample frame
        for _ in range(4096):
            item = fft_item("item")
            await self.start_item(item)
            
            item.in_real = random.randint(-32768, 32767)
            item.in_imag = random.randint(-32768, 32767)
            item.rst_n = 1
            item.valid_in = 1
            await self.finish_item(item)

        # 2.  Drop valid_in to 0 so the pipeline can drain!
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)