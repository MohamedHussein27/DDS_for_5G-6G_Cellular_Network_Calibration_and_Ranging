from pyuvm import uvm_sequence
from fft_seq_item import fft_item
import random
from pyuvm import *


class SeqNegDc(uvm_sequence):
    """
    Test Case: Maximum Negative DC Test
    Description: Pushes the absolute floor of a 16-bit signed integer (-32768) 
    continuously to test the RTL's two's-complement arithmetic, sign-extension, 
    and negative accumulation limits.
    """
    async def body(self):
        for _ in range(4096):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            # -32768 is the absolute lowest possible 16-bit signed integer
            item.in_real = -32768 
            item.in_imag = 0 
            await self.finish_item(item)
            
        # Drain the pipeline
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)

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

class fft_ten_sample_seq(uvm_sequence):
    async def body(self):
        # 1. Reset
        reset_item = fft_item("reset_item")
        reset_item.rst_n = 0
        reset_item.valid_in = 0
        await self.start_item(reset_item)
        await self.finish_item(reset_item)

        # 2. Your 9 "Interesting" Samples
        # You can change these to any values you want to test
        my_data = [100, 200, 300, 400, 500, 400, 300, 200, 100]
        
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
        for i in range(4096 - 9):
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
            # drop valid_in to 0 so the pipeline can drain safely!
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)


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
            item.in_real = random.randint(-3, 3)
            item.in_imag = random.randint(-3, 3)
            
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



class SeqImpulse(uvm_sequence):
    async def body(self):
        for i in range(4096):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = 120 if i == 0 else 0 
            item.in_imag = 0
            await self.finish_item(item)
            
        # --- ADD THIS PIPELINE DRAIN! ---
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)