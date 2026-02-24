import asyncio
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)

log = logging.getLogger(__name__)

class Consumer:
    def __init__(self):
        # Queues will be connected later
        self.in_queue = asyncio.Queue()
        self.out_queue = asyncio.Queue()
        self.start_event = None

    async def run(self):
        log.info("[Consumer] Waiting for start event...")

        # Wait for synchronization event
        await self.start_event.wait()

        log.info("[Consumer] Start event received!")

        while True:
            trx = await self.in_queue.get()
            
            await asyncio.sleep(2)    
            trx.data = trx.data * 2
            
            log.info(f"[Consumer] Processed tid={trx.tid}. Forwarding...")
            await self.out_queue.put(trx)
