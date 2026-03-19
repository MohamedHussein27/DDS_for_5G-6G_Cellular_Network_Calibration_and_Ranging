import pyuvm
from pyuvm import *
from alu_sequence_item import alu_sequence_item
from alu_pkg import opcode_t

class rst_sequence(uvm_sequence):
    async def body(self):
        # 5 clean reset pulses is plenty to verify reset logic
        for _ in range(5):
            req = alu_sequence_item.create("req")
            await self.start_item(req)
            
            req.a = 0
            req.b = 0
            req.op = opcode_t.ADD_op
            
            # Force reset to active (0)
            req.rst_n = 0 
            
            await self.finish_item(req)