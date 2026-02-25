import asyncio
import time

from producer import producer
from consumer import consumer
from checker import checker
from event_trigger import event_trigger


async def main():

    producer_to_consumer = asyncio.Queue()
    consumer_to_checker = asyncio.Queue()

    start_event = asyncio.Event()

    tasks = [
        asyncio.create_task(producer(producer_to_consumer)),
        asyncio.create_task(consumer(producer_to_consumer,
                            consumer_to_checker, start_event)),
        asyncio.create_task(checker(consumer_to_checker)),
        asyncio.create_task(event_trigger(start_event))
    ]

    await tasks[0]  # wait producer

    await producer_to_consumer.join()
    await consumer_to_checker.join()

    for t in tasks:
        t.cancel()

    print(f"[{time.strftime('%X')}] Simulation finished cleanly.")


if __name__ == "__main__":
    asyncio.run(main())
