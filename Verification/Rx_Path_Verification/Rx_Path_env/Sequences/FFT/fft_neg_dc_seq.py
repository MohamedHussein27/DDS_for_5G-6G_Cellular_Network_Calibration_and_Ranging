from pyuvm import uvm_sequence
from fft_seq_item import fft_item

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