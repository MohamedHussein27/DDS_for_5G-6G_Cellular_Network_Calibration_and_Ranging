import asyncio
import time


async def checker(out_queue):

    processed_count = 0

    while processed_count < 10:

        try:
            tx = await asyncio.wait_for(out_queue.get(), timeout=3.0)

            expected = (tx.tid * 5) * 2

            if tx.data == expected:
                print(
                    f"[{time.strftime('%X')}] Checker: PASS - TX {tx.tid} (Received: {tx.data})")
            else:
                print(f"[{time.strftime('%X')}] Checker: FAIL - TX {tx.tid}")

            out_queue.task_done()
            processed_count += 1

        except asyncio.TimeoutError:
            print(
                f"[{time.strftime('%X')}] Checker: WARNING - Timeout waiting for transaction")
