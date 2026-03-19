import pyuvm
from pyuvm import *
from alu_sequence_item import alu_sequence_item
from alu_pkg import opcode_t

class add_sequence(uvm_sequence):
    async def body(self):
        # Loop 500 times to check all ADD combinations
        for _ in range(500):
            req = alu_sequence_item.create("req")
            await self.start_item(req)
            
            # 1. Randomize A, B, and rst_n
            req.randomize() 
            
            # 2. Force the operation to ADD (overriding the randomizer)
            req.op = opcode_t.ADD_op 
            
            await self.finish_item(req)