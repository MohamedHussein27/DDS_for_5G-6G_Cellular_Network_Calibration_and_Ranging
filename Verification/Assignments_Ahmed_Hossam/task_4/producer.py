# producer.py
import asyncio
import time
from transaction import Transaction

async def producer(queue, start_time):
    """Generates transactions periodically."""
    for i in range(1, 11):
        transaction = Transaction(tid=i, data=i)
        await queue.put(transaction)
        
        current_time = time.time() - start_time
        print(f"[{current_time:.2f}s] Producer: Pushed Transaction(tid={transaction.tid}, data={transaction.data})")
        
        # Fixed 1-second delay
        await asyncio.sleep(1)