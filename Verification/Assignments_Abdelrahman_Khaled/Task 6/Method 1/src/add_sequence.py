import cocotb
from transaction import ALUTransaction


class AddSequence:
    """Sequence that generates ADD operations (op = 0)."""

    async def run(self, sequencer, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction()

            # Constrained randomization for ADD
            tr.randomize_with(lambda t: t.op == 0)

            await sequencer.send(tr)
