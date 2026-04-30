from pyuvm import *
from ifft_item import ifft_item
from ifft_seq_utils import N, _f2q

class dc_sequence(uvm_sequence):
    def __init__(self, name="dc_sequence", multi_frame=False):
        super().__init__(name)
        self.multi_frame = multi_frame

    async def body(self):
        re_val = _f2q(1.0)
        im_val = 0

        # ── First N cycles: drive the input frame ──
        for _ in range(N):
            self.seq_item = ifft_item.create("seq_item")
            await self.start_item(self.seq_item)
            
            self.seq_item.randomize_with(
                lambda rst_n, valid_in: rst_n == 1 and valid_in == 1
            )
            self.seq_item.in_real = re_val
            self.seq_item.in_imag = im_val
            
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