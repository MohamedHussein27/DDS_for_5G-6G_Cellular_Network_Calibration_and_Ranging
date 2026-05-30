from pyuvm import *
from ifft_item import ifft_item

class reset_sequence(uvm_sequence):
    def __init__(self, name="reset_sequence", multi_frame=False):
        super().__init__(name)
        self.multi_frame = multi_frame

    async def body(self):
        self.seq_item = ifft_item.create("seq_item")
        await self.start_item(self.seq_item)

        self.seq_item.randomize_with(
            lambda rst_n, valid_in: 
            rst_n == 0 and valid_in == 0
        )
        await self.finish_item(self.seq_item)