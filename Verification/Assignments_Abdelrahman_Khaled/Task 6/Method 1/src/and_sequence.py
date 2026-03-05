import cocotb
from transaction import ALUTransaction


class AndSequence:
    """Sequence that generates AND operations (op = 2)."""

    async def run(self, sequencer, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction()

            # Constrained randomization for AND
            tr.randomize_with(lambda t: t.op == 2)

            await sequencer.send(tr)
