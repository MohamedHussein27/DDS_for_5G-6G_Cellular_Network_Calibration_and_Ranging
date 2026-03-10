import random
from transaction import Transaction
from cocotb.triggers import Event
from cocotb.queue import Queue
import cocotb

class Generator:

    def __init__(self):
        self.gen2drv = Queue()
        # Initialize gen_ack as an Event so await self.gen_ack.wait() works
        self.gen_ack = Event() 
        self.num_transactions = 20

        # to end simulation
        self.finish = Event()

    async def run(self, test_type="random_test"):
        
        for i in range(self.num_transactions):
            
            # make reset sequence in every test
            if i == 0:
                tr = self.reset_sequence()
                
            # Handle standard test sequences for all other transactions
            else:
                if test_type == "add_xor_test":
                    tr = self.add_xor_sequence()
                elif test_type == "and_or_test":
                    tr = self.and_or_sequence()
                else:
                    tr = self.random_sequence()

            # Send to driver and wait for acknowledgment
            await self.gen2drv.put(tr)
            await self.gen_ack.wait()
            self.gen_ack.clear() # Clear the event so it can be used for the next loop!

        # finishing
        self.finish.set()


    # --------------------------------------------------------
    # Sequences (Using keyword arguments to prevent OverflowErrors)
    # --------------------------------------------------------
    def add_xor_sequence(self):
        a = random.randint(0, 15)
        b = random.randint(0, 15)
        op = random.choice([0, 1])  # 0=ADD 1=XOR
        
        return Transaction(rst_n=1, a=a, b=b, op=op)

    def and_or_sequence(self):
        a = random.randint(0, 15)
        b = random.randint(0, 15)
        op = random.choice([2, 3])  # 2=AND 3=OR
        return Transaction(rst_n=1, a=a, b=b, op=op)

    def random_sequence(self):
        a = random.randint(0, 15)
        b = random.randint(0, 15)
        op = random.randint(0, 3)
        return Transaction(rst_n=1, a=a, b=b, op=op)
    
    def reset_sequence(self):
        return Transaction(rst_n=0, a=0, b=0, op=0)