import cocotb 
from cocotb.triggers import * 
from pyuvm import * 
import logging

from sequence_item  import * 

class reset_sequence(uvm_sequence):
    def __init__(self, name="uvm_sequence"):
        super().__init__(name)

    async def body(self):
        self.seq_item = Sequence_Item.create("seq_item")
        await self.start_item(self.seq_item)


        self.seq_item.randomize_with(lambda rst_n: rst_n in [0])

        await self.finish_item(self.seq_item)


class add_xor_sequence(uvm_sequence):
    def __init__(self, name="uvm_sequence"):
        super().__init__(name)
    
    async def body(self):

        for _ in range(10):
            self.seq_item = Sequence_Item.create("seq_item")
            await self.start_item(self.seq_item)

            self.seq_item.randomize_with(lambda op, rst_n: op in [0, 1] and rst_n == 1)

            await self.finish_item(self.seq_item)


class and_or_sequence(uvm_sequence):
    def __init__(self, name="uvm_sequence"):
        super().__init__(name)
    
    async def body(self):
        for _ in range(10):
            self.seq_item = Sequence_Item.create("seq_item")
            await self.start_item(self.seq_item)

            self.seq_item.randomize_with(lambda op, rst_n: op in [2, 3] and rst_n == 1)

            await self.finish_item(self.seq_item)
        

class base_sequence(uvm_sequence):
    def __init__(self, name="uvm_sequence"):
        super().__init__(name)
    
    async def body(self):
        for _ in range(20):
            self.seq_item = Sequence_Item.create("seq_item")
            await self.start_item(self.seq_item)

            self.seq_item.randomize_with(lambda rst_n: rst_n == 1)

            await self.finish_item(self.seq_item)