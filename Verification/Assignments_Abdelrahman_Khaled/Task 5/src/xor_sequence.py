import cocotb
import random
from transaction import ALUTransaction


class XorSequence:
    """Class defining the Xor operation sequence."""

    async def run(self, sequencer, num_transactions=20):
        """Generates transactions and sends them to the sequencer."""

        for _ in range(num_transactions):
            a = random.randint(0, 15)
            b = random.randint(0, 15)

            tr = ALUTransaction(a, b, 1)  # op = 1 for Xor

            # Send transaction to the Sequencer
            await sequencer.send(tr)
