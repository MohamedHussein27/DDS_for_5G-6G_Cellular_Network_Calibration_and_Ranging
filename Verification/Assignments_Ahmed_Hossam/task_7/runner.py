import os
import sys
import argparse
from pathlib import Path
from cocotb.runner import get_runner
from test import run


def test_alu():
    # 1. Get the current directory path
    proj_path = os.path.dirname(os.path.abspath(__file__))
    
    # 2. Define your hardware source files
    verilog_sources = [os.path.join(proj_path, "DUT.sv")]
    
    # 3. Choose your test name from the command line or default to 'test_all'
    # This replaces the +UVM_TESTNAME logic
    uvm_test = os.getenv("TEST_NAME", "test_all")

    run(
        verilog_sources=verilog_sources,
        toplevel="ALU",            # Your Verilog module name
        module="test",             # Your Python file (test.py)
        simulator="questa",        # Use Questa!
        sim_args=[f"+UVM_TESTNAME={uvm_test}"],
        waves=True                 # This will generate a WLF file for Questa
    )

if __name__ == "__main__":
    test()