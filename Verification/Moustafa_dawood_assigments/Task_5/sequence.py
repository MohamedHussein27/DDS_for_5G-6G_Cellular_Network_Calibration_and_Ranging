from transaction import ALUTransaction
class ALUSequence:

    def __init__(self, queue):
        self.queue = queue

    async def run(self):

        for op in range(4):
            for a in range(16):
                for b in range(16):

                    item = ALUTransaction(a, b, op)

                    await self.queue.put(item)
                    
