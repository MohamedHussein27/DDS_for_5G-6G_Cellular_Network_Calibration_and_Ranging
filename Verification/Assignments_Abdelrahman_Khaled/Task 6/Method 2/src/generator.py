import cocotb
import random
from transaction import ALUTransaction


class ALUGenerator:
    def __init__(self, gen_queue):
        self.gen_queue = gen_queue

    async def generate_reset(self, num_transactions=3):
        for _ in range(num_transactions):
            tr = ALUTransaction(random.randint(
                0, 15), random.randint(0, 15), 0, 0)
            await self.gen_queue.put(tr)

    async def generate_add(self, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction(random.randint(
                0, 15), random.randint(0, 15), 0)
            await self.gen_queue.put(tr)

    async def generate_xor(self, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction(random.randint(
                0, 15), random.randint(0, 15), 1)
            await self.gen_queue.put(tr)

    async def generate_and(self, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction(random.randint(
                0, 15), random.randint(0, 15), 2)
            await self.gen_queue.put(tr)

    async def generate_or(self, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction(random.randint(
                0, 15), random.randint(0, 15), 3)
            await self.gen_queue.put(tr)

    async def generate_random_with_reset(self, num_transactions=20):
        # LAMBDA EXAMPLE: Drive reset (0) for the first 2 cycles, then high (1)
        def get_reset(idx): return 0 if idx < 2 else 1

        for i in range(num_transactions):
            tr = ALUTransaction(
                a=random.randint(0, 15),
                b=random.randint(0, 15),
                op=random.randint(0, 3),
                rst_n=get_reset(i)
            )
            await self.gen_queue.put(tr)
