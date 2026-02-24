import asyncio
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)

log = logging.getLogger(__name__)

class Checker:
    def __init__(self):
        # Interconnects will be connected later
        self.score_queue = asyncio.Queue()
        self.done_event = None
        self.passed = 0
        self.failed = 0
        self.begin_ckecking = 0

    async def run(self):
        log.info("[Checker] Active and waiting...")

        received = 0

        while received < 10:
            try:
                await self.begin_ckecking.wait()
                # Wait for transaction with 3-second timeout
                trx = await asyncio.wait_for(
                    self.score_queue.get(),
                    timeout=3
                )

                expected = (received + 1) * 20   # <-- Correct expectation

                if trx.data == expected:
                    log.info(
                        f"[Checker] PASS: tid={trx.tid} | Data={trx.data}"
                    )
                    self.passed += 1
                else:
                    log.error(
                        f"[Checker] FAIL: tid={trx.tid} "
                        f"| Expected={expected}, Got={trx.data}"
                    )
                    self.failed += 1

                received += 1

            except asyncio.TimeoutError:
                log.warning(
                    "[Checker] WARNING: No transaction received within 3 seconds!"
                )
                # Continue waiting without stopping simulation

        self.done_event.set()

    def report_test_cases(self):
        log.info("=================================")
        log.info("      SCOREBOARD REPORT          ")
        log.info("=================================")
        log.info(f" Passed: {self.passed}")
        log.info(f" Failed: {self.failed}")
        log.info("=================================")