import pyuvm
from pyuvm import *
from seq_item import *
from driver import *
from sequencer import sequencer
from monitor import monitor
from cocotb.triggers import Timer
from cocotb.clock import Clock
from cocotb_coverage.coverage import CoverPoint, coverage_db
@CoverPoint("top.out", vname="out", bins=list(range(0, 16))) 
@CoverPoint("top.c",   vname="c",   bins=list(range(0, 2)))
@CoverPoint("top.op",  vname="op",  bins=list(range(0, 4)))
def sample(out, c, op):
    pass
class subscriber(uvm_component):
    def build_phase(self):     
        self.sub_export = uvm_analysis_export("sub_export", self)
        self.sub_fifo = uvm_tlm_analysis_fifo("sub_fifo", self)
        
    def connect_phase(self):
        self.sub_export.connect(self.sub_fifo.analysis_export)
    async def run_phase(self):
        while True:
            item = await self.sub_fifo.get()
            sample(item.out, item.c, item.op)    
        
    
    def coverage_report(self):
        # Fetch the data using the exact names from the @CoverPoint decorators
        out_cov = coverage_db["top.out"].coverage          
        out_p   = coverage_db["top.out"].cover_percentage  
        
        cocotb.log.info(f"The 'out' coverage is : {out_cov}")
        cocotb.log.info(f"The 'out' coverage percentage is : {out_p}%")
        
        coverage_db.export_to_xml("alu_coverage.xml")