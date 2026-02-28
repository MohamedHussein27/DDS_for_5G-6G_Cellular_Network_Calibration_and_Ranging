import cocotb
import random
from transaction import ALUTransaction


class OrSequence:
    """Class defining the Or operation sequence."""

    async def run(self, sequencer, num_transactions=20):
        """Generates transactions Or sends them to the sequencer."""

        for _ in range(num_transactions):
            a = random.randint(0, 15)
            b = random.randint(0, 15)

            tr = ALUTransaction(a, b, 3)  # op = 3 for Or

            # Send transaction to the Sequencer
            await sequencer.send(tr)
