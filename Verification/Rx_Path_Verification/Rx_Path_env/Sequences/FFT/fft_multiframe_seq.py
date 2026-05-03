import random
from pyuvm import *
from fft_seq_item import fft_item

class fft_two_frame_seq(uvm_sequence):
    async def body(self):
        # 2 frames * 4096 samples = 8192 total continuous samples
        TOTAL_SAMPLES = 8192
        
        for i in range(TOTAL_SAMPLES):
            item = fft_item(f"rand_sample_{i}")
            
            # Setup control signals
            item.rst_n = 1
            item.valid_in = 1
            
            # Generate random 16-bit signed integers (-32768 to 32767)
            # This perfectly matches your RTL WL=16 parameter
            item.in_real = random.randint(-32768, 32767)
            item.in_imag = random.randint(-32768, 32767)
            
            await self.start_item(item)
            await self.finish_item(item)
            
        # Optional: Wait a few cycles at the very end to let valid_in drop
        for _ in range(10):
            idle_item = fft_item("idle")
            idle_item.rst_n = 1
            idle_item.valid_in = 0
            idle_item.in_real = 0
            idle_item.in_imag = 0
            await self.start_item(idle_item)
            await self.finish_item(idle_item)