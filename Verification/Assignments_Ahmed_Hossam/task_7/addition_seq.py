import pyuvm
from pyuvm import *
from seq_item import *
from cocotb.triggers import Timer
from cocotb.clock import Clock

class addition_seq(uvm_sequence):
    async def body(self):
        for a in range(16):
            for b in range(16):
                
                req = seq_item("req") 
                await self.start_item(req)
                req.reset = 1   
                req.a = a
                req.b = b
                req.op = 0       
                await self.finish_item(req)