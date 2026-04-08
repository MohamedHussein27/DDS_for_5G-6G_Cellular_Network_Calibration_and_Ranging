import pyuvm
from pyuvm import *
from seq_item import *
from cocotb.triggers import Timer
from cocotb.clock import Clock
from cocotb_coverage.coverage import CoverPoint, coverage_db


@CoverPoint("top.out", vname="out", bins=list(range(0, 16))) 
@CoverPoint("top.c",   vname="c",   bins=list(range(0, 2)))
@CoverPoint("top.op",  vname="op",  bins=list(range(0, 4)))
def sample(out, c, op):
    pass

class subscriber(uvm_subscriber):
    def build_phase(self):     
     
        self.sub_fifo = uvm_tlm_analysis_fifo("sub_fifo", self)
        
   

    async def run_phase(self):
        while True:
            item = await self.sub_fifo.get()
            sample(item.out, item.c, item.op)    
        
    
    def report_phase(self):
        out_cov = coverage_db["top.out"].coverage          
        out_p   = coverage_db["top.out"].cover_percentage  
        
        cocotb.log.info("========================================")
        cocotb.log.info("           COVERAGE REPORT              ")
        cocotb.log.info("========================================")
        cocotb.log.info(f" 'out' Hits: {out_cov}")
        cocotb.log.info(f" 'out' Percentage: {out_p}%")
        cocotb.log.info("========================================")
        
        coverage_db.export_to_xml("alu_coverage.xml")