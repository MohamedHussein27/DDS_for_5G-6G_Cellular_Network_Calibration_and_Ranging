from pyuvm import *
from ALU_seq_item import ALU_seq_item


class ALU_random_sequence(uvm_sequence):
    async def body(self):
        # Let's generate 50 random transactions
        for _ in range(50):
            item = ALU_seq_item("item")
            await self.start_item(item)
            item.randomize()
            await self.finish_item(item)
