import asyncio
import time


async def consumer(in_queue, out_queue, start_event):

    print(f"[{time.strftime('%X')}] Consumer: Waiting for start_event...")
    await start_event.wait()

    print(f"[{time.strftime('%X')}] Consumer: Started processing.")

    while True:
        tx = await in_queue.get()

        await asyncio.sleep(2)  # Processing delay

        tx.data *= 2

        await out_queue.put(tx)

        in_queue.task_done()
