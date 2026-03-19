from pyuvm import *
from alu_sequence_item import alu_sequence_item

class random_op(uvm_sequence):
    async def body(self):
        for _ in range(20000):
            item = alu_sequence_item("item")
            await self.start_item(item)
            item.randomize()
            await self.finish_item(item)