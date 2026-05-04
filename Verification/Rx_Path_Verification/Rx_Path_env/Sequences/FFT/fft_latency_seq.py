import random
from pyuvm import *
from fft_seq_item import fft_item

class SeqLatency(uvm_sequence):
    async def body(self):
        #  OPEN A TEXT FILE FOR MATLAB
        # "w" means write mode (it will overwrite the file each time you run)
        with open("matlab_inputs.txt", "w") as f:
            
            for _ in range(4096):
                # 1. Generate the random numbers FIRST
                r_val = random.randint(-7, 7)
                i_val = random.randint(-7, 7)
                
                # 2. Write them to the text file (Format: "Real Imag\n")
                f.write(f"{r_val} {i_val}\n")
                
                # 3. Send those exact same numbers to the RTL
                item = fft_item("item")
                await self.start_item(item)
                item.in_real = r_val
                item.in_imag = i_val
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