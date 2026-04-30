import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from dds_seq_item  import * 
class dds_reset_recovery_seq(uvm_sequence):
    def __init__(self, name="dds_reset_recovery_seq"):
        super().__init__(name)

    async def body(self):
        for _ in range(5):
            # Step 1: Send a reset state to kill whatever the hardware is doing
            req_kill = dds_seq_item("req_kill")
            await self.start_item(req_kill)
            req_kill.randomize()
            req_kill.rst_n = 0  # FORCE RESET ACTIVE
            req_kill.enable = 1  # Ensure enable is high for valid transactions 
            await self.finish_item(req_kill)
            
            # Step 2: Immediately send a valid transaction to ensure it woke up
            req_recover = dds_seq_item("req_recover")
            await self.start_item(req_recover)
            req_recover.randomize()
            req_recover.enable = 1  # Ensure enable is high for valid transactions
            req_recover.rst_n = 1 # FORCE RESET INACTIVE
            await self.finish_item(req_recover)