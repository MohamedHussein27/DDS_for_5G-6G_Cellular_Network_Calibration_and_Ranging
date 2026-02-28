import cocotb
import random
from transaction import ALUTransaction


class AndSequence:
    """Class defining the And operation sequence."""

    async def run(self, sequencer, num_transactions=20):
        """Generates transactions and sends them to the sequencer."""

        for _ in range(num_transactions):
            a = random.randint(0, 15)
            b = random.randint(0, 15)

            tr = ALUTransaction(a, b, 2)  # op = 2 for And

            # Send transaction to the Sequencer
            await sequencer.send(tr)
