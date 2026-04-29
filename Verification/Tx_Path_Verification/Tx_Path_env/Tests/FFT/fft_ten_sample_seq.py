from pyuvm import *
from fft_seq_item import fft_item

class fft_ten_sample_seq(uvm_sequence):
    async def body(self):
        # 1. Reset
        reset_item = fft_item("reset_item")
        reset_item.rst_n = 0
        reset_item.valid_in = 0
        await self.start_item(reset_item)
        await self.finish_item(reset_item)

        # 2. Your 10 "Interesting" Samples
        # You can change these to any values you want to test
        my_data = [100, 200, 300, 400, 500, 400, 300, 200, 100, 50]
        
        for val in my_data:
            item = fft_item("data_item")
            item.rst_n    = 1
            item.valid_in = 1
            item.in_real  = val
            item.in_imag  = 0
            await self.start_item(item)
            await self.finish_item(item)

        # 3. ZERO-PADDING (Crucial)
        # We MUST send the remaining 4086 samples as zeros to finish the frame
        for i in range(4096 - 10):
            zero_item = fft_item(f"pad_{i}")
            zero_item.rst_n    = 1
            zero_item.valid_in = 1
            zero_item.in_real  = 0
            zero_item.in_imag  = 0
            await self.start_item(zero_item)
            await self.finish_item(zero_item)

        # 4. Wait for the pipeline to drain
        for _ in range(20):
            idle = fft_item("idle")
            idle.rst_n = 1
            idle.valid_in = 0
            await self.start_item(idle)
            await self.finish_item(idle)