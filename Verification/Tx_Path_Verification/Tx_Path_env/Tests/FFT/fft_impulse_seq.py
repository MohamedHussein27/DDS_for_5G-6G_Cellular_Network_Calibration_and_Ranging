from pyuvm import uvm_sequence
from fft_seq_item import fft_item

class SeqImpulse(uvm_sequence):
    async def body(self):
        for i in range(4096):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            
            # Blast a massive impulse so it survives the 12 stages of scaling!
            # 16384 is exactly 2^14. After dividing by 4096, it will output a clean 4.
            item.in_real = 120 if i == 0 else 0 
            item.in_imag = 0
            
            await self.finish_item(item)