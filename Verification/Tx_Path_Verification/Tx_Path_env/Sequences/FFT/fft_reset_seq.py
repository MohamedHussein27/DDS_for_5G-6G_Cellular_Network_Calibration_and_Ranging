from pyuvm import uvm_sequence
from fft_seq_item import fft_item

class FftResetSeq(uvm_sequence):
    def __init__(self, name="FftResetSeq"):
        super().__init__(name)

    async def body(self):
        # 1. Assert Reset (Keep valid_in LOW so the scoreboard ignores it)
        item = fft_item("reset_assert")
        await self.start_item(item)
        item.rst_n = 0
        item.valid_in = 0
        item.in_real = 0
        item.in_imag = 0
        await self.finish_item(item)
        
        # 2. Release Reset (Hardware is now awake and ready)
        item = fft_item("reset_release")
        await self.start_item(item)
        item.rst_n = 1
        item.valid_in = 0
        item.in_real = 0
        item.in_imag = 0
        await self.finish_item(item)