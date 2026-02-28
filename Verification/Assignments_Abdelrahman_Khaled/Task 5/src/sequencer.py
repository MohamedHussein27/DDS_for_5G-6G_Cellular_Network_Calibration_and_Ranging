import cocotb
from cocotb.queue import Queue


class ALUSequencer:
    def __init__(self):
        # Cocotb's queue, perfectly synced with the simulator time
        self.queue = Queue()

    async def send(self, transaction):
        """Receives transactions from sequences and queues them for the driver."""
        self.queue.put_nowait(transaction)
