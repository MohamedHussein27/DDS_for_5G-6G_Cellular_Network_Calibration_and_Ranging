# checker.py
import asyncio
import time

async def checker(queue, start_time):
    """Verifies correctness of processed data."""
    while True:
        try:
            # Timeout detection: 3 seconds
            transaction = await asyncio.wait_for(queue.get(), timeout=3.0)
            
            expected = transaction.data * 2
            current_time = time.time() - start_time
            
            if transaction.processed_data == expected:
                print(f"[{current_time:.2f}s] Checker: PASS for tid={transaction.tid} | Expected: {expected}, Got: {transaction.processed_data}")
            else:
                print(f"[{current_time:.2f}s] Checker: FAIL for tid={transaction.tid}")
                
            queue.task_done()
            
        except asyncio.TimeoutError:
            print(f"[{time.time() - start_time:.2f}s] Checker: WARNING - No transaction received within 3 seconds!")