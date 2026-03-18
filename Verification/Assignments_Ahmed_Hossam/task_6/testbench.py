import os
import cocotb
from environment import *
from cocotb.triggers import * 
from cocotb.clock    import Clock
from cocotb_coverage.coverage import coverage_db

@cocotb.test()
async def tb(dut):
    # Read the variable here as well for logging
    seq_name = os.environ.get("TEST_SEQ", "all")
    
    cocotb.log.info(f" === STARTING SIMULATION | TEST_SEQ: {seq_name} === ")
    
    e = environment()
    clk = Clock(dut.clk, 10, units="ns")
    await cocotb.start(clk.start())
    
    cocotb.start_soon(e.run_environment(dut))
    
    await e.join_any.wait()
    coverage_db.export_to_xml(filename="ALU_coverage.xml")