from pyuvm import *
from ALU_seq_item import ALU_seq_item


class ALU_and_sequence(uvm_sequence):
    async def body(self):
        for a in range(16):
            for b in range(16):
                item = ALU_seq_item("item")
                await self.start_item(item)
                item.op = 2  # 2'b10 [cite: 34]
                item.a = a
                item.b = b
                await self.finish_item(item)
