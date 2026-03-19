import pyuvm
from pyuvm import *
from seq_item import *
from cocotb.triggers import Timer
from cocotb.clock import Clock

class driver(uvm_driver):   
    async def run_phase(self):
        while True: 
            
            seq_item = await self.seq_item_port.get_next_item()
            
            
            self.dut.reset.value = seq_item.reset
            self.dut.a.value = seq_item.a
            self.dut.b.value = seq_item.b
            self.dut.op.value = seq_item.op
            
            
            await Timer(10, units="ns")   
            
            self.seq_item_port.item_done()