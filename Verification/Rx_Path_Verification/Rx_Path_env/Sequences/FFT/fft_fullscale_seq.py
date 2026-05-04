from pyuvm import uvm_sequence
from fft_seq_item import fft_item

class SeqFullScale(uvm_sequence):
    async def body(self):
        # 1. Blast the alternating Nyquist maximums
        for i in range(4096):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = 32767 if i % 2 == 0 else -32768
            item.in_imag = 32767 if i % 2 == 0 else -32768
            await self.finish_item(item)
            
        # 2.  Drop valid_in to 0 so the pipeline drains safely!
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)