from cocotb.triggers import with_timeout
from cocotb.result import SimTimeoutError

class Checker:
    def __init__(self, in_queue, out_queue):
        self.in_queue = in_queue
        self.out_queue = out_queue

    async def run(self):
        while True:
            try:
                
                txn = await with_timeout(self.in_queue.get(), 30, 'ns')
                result = await with_timeout(self.out_queue.get(), 30, 'ns')

                expected = txn.data * 2

                if result == expected:
                    print(f"[CHECKER] PASS TID={txn.tid}")
                else:
                    print(f"[CHECKER] FAIL TID={txn.tid}")

            except SimTimeoutError: 
                print("[CHECKER WARNING] Timeout detected")