import logging
import asyncio
from environment import Environment


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)

log = logging.getLogger(__name__)

async def main():
    # 1. Build Phase
    env = Environment()
    log.info("Simulation begins")
    
    # 2. Connect Phase
    env.connect()
    
    # 3. Run Phase
    await env.run_environment()
    
    log.info("Waiting for the Checker to signal completion...")
    
    # Keep the simulation alive until the join_any event is triggered
    await env.join_any.wait()
    
    # 4. Report Phase
    log.info("Completion signal received. Ending test.")
    env.s.report_test_cases()

if __name__ == "__main__":
    asyncio.run(main())