from transaction import ALUTransaction

class BaseSequence:
    def __init__(self, count=256): # 16 x 16 = 256 items per operation
        self.seq2drv = None
        self.count = count

    def connect(self, mailbox):
        self.seq2drv = mailbox

class AddSequence(BaseSequence):
    async def run(self):
        for a in range(16):
            for b in range(16):
                item = ALUTransaction(a, b, 0)
                await self.seq2drv.put(item)

class XorSequence(BaseSequence):
    async def run(self):
        for a in range(16):
            for b in range(16):
                item = ALUTransaction(a, b, 1)
                await self.seq2drv.put(item)

class AndSequence(BaseSequence):
    async def run(self):
        for a in range(16):
            for b in range(16):
                item = ALUTransaction(a, b, 2)
                await self.seq2drv.put(item)

class OrSequence(BaseSequence):
    async def run(self):
        for a in range(16):
            for b in range(16):
                item = ALUTransaction(a, b, 3)
                await self.seq2drv.put(item)
