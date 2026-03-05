import cocotb
from transaction import ALUTransaction


class OrSequence:
    """Sequence that generates OR operations (op = 3)."""

    async def run(self, sequencer, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction()

            # Constrained randomization for OR
            tr.randomize_with(lambda t: t.op == 3)

            await sequencer.send(tr)
