from pyuvm import uvm_sequence
from fft_seq_item import FftSeqItem

class SeqZero(uvm_sequence):
    async def body(self):
        for _ in range(4096):
            item = FftSeqItem("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = 0
            item.in_imag = 0
            await self.finish_item(item)