# main.py
import asyncio
import time

# Import the coroutines from your other files
from producer import producer
from consumer import consumer
from checker import checker

async def trigger_start(start_event, start_time):
    """Triggers the start event after 5 seconds."""
    await asyncio.sleep(5)
    print(f"[{time.time() - start_time:.2f}s] Trigger: Firing start_event now!")
    start_event.set()

async def main():
    start_time = time.time()
    
    # Create the communication channels
    in_queue = asyncio.Queue()
    out_queue = asyncio.Queue()
    start_event = asyncio.Event()

    # Launch all coroutines concurrently, passing the start_time to each
    prod_task = asyncio.create_task(producer(in_queue, start_time))
    cons_task = asyncio.create_task(consumer(in_queue, out_queue, start_event, start_time))
    check_task = asyncio.create_task(checker(out_queue, start_time))
    trigger_task = asyncio.create_task(trigger_start(start_event, start_time))

    # Wait for the simulation to finish cleanly
    await prod_task
    await in_queue.join()
    await out_queue.join()

    # Shut down the infinite loops
    cons_task.cancel()
    check_task.cancel()
    
    print(f"[{time.time() - start_time:.2f}s] Simulation Controller: Simulation finished cleanly.")

if __name__ == "__main__":
    # This is the single entry point to run your whole environment
    asyncio.run(main())