import cocotb
import random
from transaction import ALUTransaction


class RandomSequence:
    """Class defining the random operation sequence."""

    async def run(self, sequencer, num_transactions=20):
        """Generates transactions random sends them to the sequencer."""

        for _ in range(num_transactions):
            a = random.randint(0, 15)
            b = random.randint(0, 15)
            op = random.randint(0, 3)

            tr = ALUTransaction(a, b, op)

            # Send transaction to the Sequencer
            await sequencer.send(tr)
