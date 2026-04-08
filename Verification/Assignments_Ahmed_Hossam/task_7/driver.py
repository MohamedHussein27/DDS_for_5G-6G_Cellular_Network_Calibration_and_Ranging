import pyuvm
from pyuvm import *
from seq_item import *
import cocotb
from cocotb.triggers import Timer
from cocotb.clock import Clock
from cocotb.triggers import FallingEdge

class driver(uvm_driver):   
    async def run_phase(self):
        while True: 
            
            seq_item = await self.seq_item_port.get_next_item()
            
            
            cocotb.top.reset.value = seq_item.reset
            cocotb.top.a.value = seq_item.a
            cocotb.top.b.value = seq_item.b
            cocotb.top.op.value = seq_item.op
            
            
            await FallingEdge(cocotb.top.clk) 
            
            self.seq_item_port.item_done()