import pyuvm
from pyuvm import *
from seq_item import *
import cocotb
from cocotb.triggers import RisingEdge, ReadOnly 

class monitor(uvm_monitor):
    def build_phase(self):
        self.mon_port = uvm_analysis_port("mon_port", self)

    async def run_phase(self):
        while True:
            await RisingEdge(cocotb.top.clk)
            
            
            await ReadOnly() 
            
            item = seq_item()
            
            item.reset = cocotb.top.reset.value
            item.a = cocotb.top.a.value
            item.b = cocotb.top.b.value
            item.op = cocotb.top.op.value
            item.c = cocotb.top.c.value
            item.out = cocotb.top.out.value
            
            self.mon_port.write(item)