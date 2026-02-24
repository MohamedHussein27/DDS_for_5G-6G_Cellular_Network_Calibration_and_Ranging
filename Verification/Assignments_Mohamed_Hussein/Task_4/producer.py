import asyncio
from transaction import Transaction

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)

log = logging.getLogger(__name__)


class Producer:
    def __init__(self):
        # The queue will be connected later by the environment
        self.drive_q = asyncio.Queue()

    async def run(self):
        log.info("[Producer] Starting generation...")
        for i in range(1, 11):
            await asyncio.sleep(1)
            
            trx = Transaction(tid=i, data=i * 10)
            log.info(f"[Producer] Generated: tid={trx.tid}, data={trx.data}")
            
            # Put data into the dynamically connected queue
            await self.drive_q.put(trx)
        
        log.info("[Producer] Finished generating 10 transactions.")
