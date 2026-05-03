"""
    Module: custom_data_seq.py
"""
import cocotb
from pyuvm import *
from fft_seq_item import fft_item
import os
from pyuvm import uvm_sequence
from fft_seq_item import fft_item

class SeqCustomData(uvm_sequence):
    """
    Reads a 1D array of integers from a text file and streams 
    them into the real port of the FFT. The imaginary port is held at 0.
    """
    def __init__(self, name="SeqCustomData"):
        super().__init__(name)

    async def body(self):
        filename = "rtl_dds_out_testing.txt"
        
        # Check if file exists
        if not os.path.exists(filename):
            cocotb.log.error(f"Cannot find {filename}! Please ensure the file is in the simulation directory.")
            return

        # 1. Read all samples from the file
        samples = []
        with open(filename, "r") as f:
            for line in f:
                # Strip whitespace and ignore empty lines
                val_str = line.strip()
                if val_str:
                    try:
                        samples.append(int(val_str))
                    except ValueError:
                        cocotb.log.warning(f"Ignored non-integer line in {filename}: {line}")
        
        cocotb.log.info(f"Loaded {len(samples)} samples from {filename}.")
        
        # 2. Stream the samples into the FFT
        for i, val in enumerate(samples):
            item = fft_item(f"item_{i}")
            await self.start_item(item)
            
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = val  # The data from your text file
            item.in_imag = 0    # Hardcoded to 0 for a real-only signal
            
            await self.finish_item(item)
            
        # 3. ZERO-PADDING (Crucial for FFT)
        # If your text file has fewer than 4096 samples, you MUST pad the rest with zeros 
        # so the FFT completes its frame and drops the valid_out flag.
        samples_sent = len(samples)
        if samples_sent < 4096:
            pad_length = 4096 - samples_sent
            cocotb.log.info(f"Zero-padding the remaining {pad_length} samples to complete the frame...")
            for i in range(pad_length):
                pad_item = fft_item(f"pad_{i}")
                await self.start_item(pad_item)
                pad_item.rst_n = 1
                pad_item.valid_in = 1
                pad_item.in_real = 0
                pad_item.in_imag = 0
                await self.finish_item(pad_item)

        # 4. Drain the Pipeline
        cocotb.log.info("Frame complete. Dropping valid_in to drain pipeline.")
        idle_item = fft_item("idle")
        await self.start_item(idle_item)
        idle_item.rst_n = 1
        idle_item.valid_in = 0
        idle_item.in_real = 0
        idle_item.in_imag = 0
        await self.finish_item(idle_item)