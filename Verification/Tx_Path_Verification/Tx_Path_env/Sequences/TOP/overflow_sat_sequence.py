import pyuvm
from pyuvm import *
from top_seq_item import top_item
import cocotb
import random
class overflow_saturation_seq(uvm_sequence):
    """
    TC-007: Overflow / Saturation Behaviour
    Stimulus: Drive worst-case maximum-amplitude stimulus to force the 
    internal adders/mixers to hit the +/-32767 rails.
    """
    def __init__(self, name="overflow_saturation_seq"):
        super().__init__(name)
        
    async def body(self):
        req = top_item.create("req")
        await self.start_item(req)
        
        req.randomize()
        
        req.rst_n = 1
        req.dds_enable = 1
        req.ofdm_rd_en = 1
        req.cycles = 4096 
        
        # Force worst-case maximum amplitudes on the OFDM inputs.
        # Combined with the DDS output, this is mathematically guaranteed 
        # to cause an internal bit-growth overflow, testing the clipping logic.
        extremes = [-32768, 32767]
        req.ofdm_in_real = random.choice(extremes)
        req.ofdm_in_imag = random.choice(extremes)
        
        await self.finish_item(req)