import os
import sys
import argparse
from pathlib import Path
from cocotb.runner import get_runner


def get_questa_runner(testcase=None):
    SIM = os.getenv("SIM", "questa")
    project_path = Path(__file__).parent.resolve()

    # Make sure Python can find your testbench
    sys.path.insert(0, str(project_path / "src"))

    # Path to your SystemVerilog design
    alu_file = project_path / "src" / "dut.sv"

    # Get the simulator runner
    runner = get_runner(SIM)

    # Build step — compile your DUT
    runner.build(
        sources=[alu_file],
        hdl_toplevel="ALU",
        always=True,
        clean=True,
        waves=True,
    )

    # Test step — run cocotb
    runner.test(
        test_module="testbench",
        testcase=testcase,          # Pass the specific test case here
        hdl_toplevel="ALU",
        hdl_toplevel_lang="verilog",
        verbose=True,
        gui=True,                   # Set to False if you want to run purely in terminal
        test_args=[
            "-l", "transcript.log"
        ],
    )


if __name__ == "__main__":
    # Add an argument parser to select tests from the command line
    parser = argparse.ArgumentParser(description="Run Cocotb Tests")
    parser.add_argument("--test", type=str, default=None,
                        help="Name of specific test to run (e.g., run_add_test)")
    args = parser.parse_args()

    get_questa_runner(testcase=args.test)
