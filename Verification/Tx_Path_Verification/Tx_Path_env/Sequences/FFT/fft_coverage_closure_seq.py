from pyuvm import uvm_sequence
from fft_seq_item import fft_item
import logging
import cocotb
from cocotb.triggers import RisingEdge
from pyuvm import *
import pyuvm

class FftCoverageClosureSeq(uvm_sequence):
    """
    Targeted sequence to close 16-bit extreme corner cases.
    Injects all 25 permutations of MAX_VAL, MIN_VAL, ZERO, and mid-range values,
    followed by pure DC saturation frames.
    """
    async def body(self):
        WL = 16
        MAX_VAL = (1 << (WL - 1)) - 1  # 32767
        MIN_VAL = -(1 << (WL - 1))     # -32768
        N = 4096
        
        # The 5 regions defined in your coverage model
        val_map = [
            MIN_VAL,    # min_neg
            -16384,     # negative
            0,          # zero
            16383,      # positive
            MAX_VAL     # max_pos
        ]
        
        # Generate all 25 cross combinations (real, imag)
        corner_pairs = []
        for r_val in val_map:
            for i_val in val_map:
                corner_pairs.append((r_val, i_val))
        
        cocotb.log.info("Starting Coverage Closure: Injecting 25-Permutation Matrix...")
        
        # ==========================================
        # PHASE 1: Back-to-Back Permutation Frames 
        # (8192 samples total to hit valid_out crosses)
        # ==========================================
        for i in range(N * 2):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            
            # Loop endlessly through the 25 permutations
            pair = corner_pairs[i % 25]
            item.in_real = pair[0]
            item.in_imag = pair[1]
            
            await self.finish_item(item)

        # ==========================================
        # PHASE 2: Output Saturation (Max Positive DC)
        # ==========================================
        cocotb.log.info("Starting Coverage Closure: Max Positive Saturation...")
        for i in range(N):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = MAX_VAL
            item.in_imag = MAX_VAL
            await self.finish_item(item)

        # ==========================================
        # PHASE 3: Output Saturation (Min Negative DC)
        # ==========================================
        cocotb.log.info("Starting Coverage Closure: Min Negative Saturation...")
        for i in range(N):
            item = fft_item("item")
            await self.start_item(item)
            item.rst_n = 1
            item.valid_in = 1
            item.in_real = MIN_VAL
            item.in_imag = MIN_VAL
            await self.finish_item(item)
            
        # ==========================================
        # PHASE 4: PIPELINE DRAIN
        # ==========================================
        cocotb.log.info("Starting Coverage Closure: Final Pipeline Flush...")
        for i in range(N + 100):
            idle_item = fft_item("idle")
            await self.start_item(idle_item)
            idle_item.rst_n = 1
            idle_item.valid_in = 0
            idle_item.in_real = 0
            idle_item.in_imag = 0
            await self.finish_item(idle_item)