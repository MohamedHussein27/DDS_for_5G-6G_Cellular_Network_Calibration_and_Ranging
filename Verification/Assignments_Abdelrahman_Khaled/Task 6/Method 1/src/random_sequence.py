import cocotb
from transaction import ALUTransaction


class RandomSequence:
    """Sequence that generates random operations with reset high."""

    async def run(self, sequencer, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction()

            # Constrain rst_n to be 1 (No Reset) for normal random operations
            tr.randomize_with(lambda t: t.rst_n == 1)

            await sequencer.send(tr)
