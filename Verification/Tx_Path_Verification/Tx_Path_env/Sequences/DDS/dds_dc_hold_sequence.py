import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_dc_hold_seq(uvm_sequence):
    """Tests if the DDS can hold a stable DC amplitude without drifting."""
    def __init__(self, name="dds_dc_hold_seq"):
        super().__init__(name)

    async def body(self):
        req = dds_seq_item.create("req")
        req.rst_n = 1
        req.enable = 1
        req.FTW_start = 0x20000000 # 45-degree phase offset
        req.FTW_step = 0           # NO movement allowed
        req.cycles = 500           # Hold it steady for 500 cycles
        for _ in range(req.cycles+1):
            await self.start_item(req)
            await self.finish_item(req)
      
        
        # Safely power down
        req_stop = dds_seq_item.create("req_stop")
        await self.start_item(req_stop)
        req_stop.rst_n = 1
        req_stop.enable = 0
        req_stop.cycles = 1
        await self.finish_item(req_stop)      
        req_stop = dds_seq_item.create("req_stop")
        await self.start_item(req_stop)
        req_stop.rst_n = 1
        req_stop.enable = 0
        req_stop.cycles = 1
        await self.finish_item(req_stop)