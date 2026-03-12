import os  
from transaction import transaction
from cocotb.triggers import *
import cocotb
import cocotb.queue

class generator():
    def __init__(self, join_any, name="GENERATOR"): 
       self.name         = name
       self.transac      = transaction()
       self.gen_mail     = cocotb.queue.Queue()
       self.gen_handover = Event(name=None) 
       self.join_any     = join_any

    async def run_generator(self, dut_generator): 
        iteration_number = 1000
        
        # Grab the variable from the Makefile. If it's not set, default to "all"
        seq_name = os.environ.get("TEST_SEQ", "all").lower()
        
        cocotb.log.info(f"[Generator] Starting sequence generation mode: {seq_name}")

        for i in range(iteration_number): 
            self.gen_handover.clear() 
            
            # Always run reset on the first iteration
            if i == 0:
                await self.reset_sequence()
            else:
                # Route to the specific sequence requested by the Makefile
                if seq_name == "add":
                    await self.add()
                elif seq_name == "xor":
                    await self.xor()
                elif seq_name == "and":
                    await self.anding()
                elif seq_name == "or":
                    await self.oring()
                else: 
                    # The default "all" behavior (runs everything sequentially)
                    await self.add()
                    await self.xor()
                    await self.anding()
                    await self.oring()

        await FallingEdge(dut_generator.clk)
        self.join_any.set()

    """ ************* **************** Sequence Generation ****************** ***************"""
    async def reset_sequence (self):
        self.transac = transaction()
        self.transac.randomize_with(lambda reset: reset == 0  )
        cocotb.log.info("[Generator] Sending To The Driver..... ") 
        await self.gen_mail.put(self.transac)  
        await self.gen_handover.wait()      

    async def add(self) :
        self.gen_handover.clear() 
        self.transac = transaction()
        self.transac.randomize_with(lambda reset :reset == 1 ,lambda op : op ==0 )
        cocotb.log.info("[Generator] Sending To The Driver..... ")
        await self.gen_mail.put(self.transac) 
        await self.gen_handover.wait()

    async def xor(self) :
        self.gen_handover.clear() 
        self.transac = transaction()
        self.transac.randomize_with(lambda reset :reset == 1 ,lambda op :op == 1 )
        cocotb.log.info("[Generator] Sending To The Driver..... ")
        await self.gen_mail.put(self.transac) 
        await self.gen_handover.wait()
     
    async def anding(self) :
        self.gen_handover.clear() 
        self.transac = transaction()
        self.transac.randomize_with(lambda reset :reset == 1 ,lambda op :op == 2 )
        cocotb.log.info("[Generator] Sending To The Driver..... ")
        await self.gen_mail.put(self.transac) 
        await self.gen_handover.wait()
    
    async def oring(self) :
        self.gen_handover.clear() 
        self.transac = transaction()
        self.transac.randomize_with(lambda reset :reset == 1 ,lambda op :op == 3 )
        cocotb.log.info("[Generator] Sending To The Driver..... ")
        await self.gen_mail.put(self.transac) 
        await self.gen_handover.wait()