# consumer.py
import asyncio
import time

async def consumer(in_queue, out_queue, start_event, start_time):
    """Receives transactions and processes them with delay."""
    print(f"[{time.time() - start_time:.2f}s] Consumer: Waiting for start_event to begin processing...")
    await start_event.wait()
    print(f"[{time.time() - start_time:.2f}s] Consumer: start_event received! Processing started.")

    while True:
        transaction = await in_queue.get()
        
        # Simulate processing delay
        await asyncio.sleep(2) 
        
        # Modify data
        transaction.processed_data = transaction.data * 2
        
        print(f"[{time.time() - start_time:.2f}s] Consumer: Processed tid={transaction.tid}")
        
        await out_queue.put(transaction)
        in_queue.task_done()