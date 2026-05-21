from pyuvm import uvm_sequence
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

        # ==========================================
        # PIPELINE DRAIN
        # ==========================================
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)