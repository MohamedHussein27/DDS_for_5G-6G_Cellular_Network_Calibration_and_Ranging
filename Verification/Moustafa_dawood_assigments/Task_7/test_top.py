import cocotb
from cocotb.clock import Clock
from pyuvm import *
import alu_test  # <--- CRITICAL: This loads the tests into the UVM Factory

@cocotb.test()
async def test_alu(dut):
    cocotb.start_soon(Clock(dut.clk, 2, unit="ns").start())
    
    # ADD THESE TWO LINES:
    print(f">>> ALL PLUSARGS: {dict(cocotb.plusargs)}")
    test_name = cocotb.plusargs.get("UVM_TESTNAME") or "regression_test"
    print(f">>> RUNNING TEST: {test_name}")
    
    await uvm_root().run_test(test_name)