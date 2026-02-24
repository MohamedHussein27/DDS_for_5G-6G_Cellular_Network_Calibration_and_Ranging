import asyncio
from producer import Producer
from consumer import Consumer
from checker import Checker

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)

log = logging.getLogger(__name__)


class Environment:
    def __init__(self):
        # events
        self.join_any = asyncio.Event()
        self.start_event = asyncio.Event()
        
        # instances
        self.p = Producer()
        self.c = Consumer()
        self.s = Checker()

    def connect(self):
        # --- UVM-style "connect_phase" ---
        log.info("[Environment] Connecting components...")
        
        # Wire up the Consumer
        self.c.in_queue = self.p.drive_q
        
        # Wire up the Checker
        self.s.score_queue = self.c.out_queue

        # connecting events
        self.s.done_event = self.join_any
        self.c.start_event = self.start_event
        self.p.begin_driving = self.start_event
        self.s.begin_ckecking = self.start_event

    async def run_environment(self):
        # --- UVM-style "run_phase" ---
        log.info("[Environment] Launching tasks...")
        asyncio.create_task(self.p.run())
        asyncio.create_task(self.c.run())
        asyncio.create_task(self.s.run())
        asyncio.create_task(self.trigger_start()) # launching start event

    async def trigger_start(self):
        log.info("[Environment] Triggering start event in 5 seconds...")
        await asyncio.sleep(5)
        self.start_event.set()
        log.info("[Environment] Start event triggered!")