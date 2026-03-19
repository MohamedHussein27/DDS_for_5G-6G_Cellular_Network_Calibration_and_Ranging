import pyuvm
from pyuvm import *
from seq_item import *
from cocotb.triggers import Timer
from cocotb.clock import Clock
class monitor(uvm_monitor):
    def build_phase(self):
        self.mon_port = uvm_analysis_port("mon_port")
    async def run_phase(self):
        while True:
            await Timer(10, units="ns") 
            item = seq_item()
            item.reset = self.dut.reset.value
            item.a = self.dut.a.value
            item.b = self.dut.b.value
            item.op = self.dut.op.value
            item.c = self.dut.c.value
            item.out = self.dut.out.value
            
            
            self.mon_port.write(item)
            