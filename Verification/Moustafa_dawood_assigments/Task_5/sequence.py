import random
import cocotb
from transaction import ALUTransaction

# ==========================================================
# ADD SEQUENCE (op = 0)
# ==========================================================
class AddSequence:
    def __init__(self, queue, count=500):
        self.queue = queue
        self.count = count

    async def run(self):
        for _ in range(self.count):
            a = random.randint(0, 15)
            b = random.randint(0, 15)
            # 0 is the OpCode for ADD
            item = ALUTransaction(a, b, 0)
            self.queue.put_nowait(item)


# XOR SEQUENCE (op = 1)

class XorSequence:
    def __init__(self, queue, count=500):
        self.queue = queue
        self.count = count

    async def run(self):
        for _ in range(self.count):
            a = random.randint(0, 15)
            b = random.randint(0, 15)
            # 1 is the OpCode for XOR
            item = ALUTransaction(a, b, 1)
            self.queue.put_nowait(item)


# AND SEQUENCE (op = 2)

class AndSequence:
    def __init__(self, queue, count=500):
        self.queue = queue
        self.count = count

    async def run(self):
        for _ in range(self.count):
            a = random.randint(0, 15)
            b = random.randint(0, 15)
            # 2 is the OpCode for AND
            item = ALUTransaction(a, b, 2)
            self.queue.put_nowait(item)


# OR SEQUENCE (op = 3)

class OrSequence:
    def __init__(self, queue, count=500):
        self.queue = queue
        self.count = count

    async def run(self):
        for _ in range(self.count):
            a = random.randint(0, 15)
            b = random.randint(0, 15)
            # 3 is the OpCode for OR
            item = ALUTransaction(a, b, 3)
            self.queue.put_nowait(item)


