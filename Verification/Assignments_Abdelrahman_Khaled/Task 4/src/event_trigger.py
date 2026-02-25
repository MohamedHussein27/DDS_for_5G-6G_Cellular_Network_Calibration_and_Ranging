import asyncio
import time


async def event_trigger(start_event):

    await asyncio.sleep(5)

    print(f"[{time.strftime('%X')}] Event Trigger: Activating Consumer!")

    start_event.set()
