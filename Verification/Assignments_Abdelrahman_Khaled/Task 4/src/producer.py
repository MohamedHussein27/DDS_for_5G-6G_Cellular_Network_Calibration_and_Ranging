import asyncio
import time
from transaction import Transaction


async def producer(queue):

    for i in range(1, 11):  # Generate 10 transactions
        await asyncio.sleep(1)  # Fixed 1 second delay

        tx = Transaction(tid=i, data=i * 5)

        await queue.put(tx)

        print(
            f"[{time.strftime('%X')}] Producer: Sent TX {tx.tid} with data {tx.data}")
