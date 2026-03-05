import cocotb
from transaction import ALUTransaction


class XorSequence:
    """Sequence that generates XOR operations (op = 1)."""

    async def run(self, sequencer, num_transactions=20):
        for _ in range(num_transactions):
            tr = ALUTransaction()

            # Constrained randomization for XOR
            tr.randomize_with(lambda t: t.op == 1)

            await sequencer.send(tr)
