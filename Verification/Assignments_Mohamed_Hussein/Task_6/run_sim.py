import os
import sys
import argparse
from pathlib import Path
from cocotb.runner import get_runner

def get_questa_runner(testcase=None, run_regression=False, use_gui=True):
    SIM = os.getenv("SIM", "questa")
    project_path = Path(__file__).parent.resolve()

    # Make sure Python can find your testbench (test_alu.py) in the current directory
    sys.path.insert(0, str(project_path))

    # Path to your SystemVerilog design (Using .as_posix() to safely handle Windows spaces)
    alu_file = project_path.joinpath("DUT.sv").as_posix()

    # Get the simulator runner
    runner = get_runner(SIM)

    # Build step — compile your DUT
    runner.build(
        sources=[alu_file],
        hdl_toplevel="ALU",
        always=True,
        clean=True,
        build_args=["-sv"] # Pass SystemVerilog flag to Questa
    )

    # Determine which tests to run based on the arguments
    if run_regression:
        target_tests = ["add_xor_test", "and_or_test", "random_test"]
        print("\n=== RUNNING REGRESSION SUITE ===")
    elif testcase:
        target_tests = [testcase]
        print(f"\n=== RUNNING SPECIFIC TEST: {testcase} ===")
    else:
        target_tests = ["random_test"] # Default fallback
        print("\n=== RUNNING DEFAULT TEST: random_test ===")

    # Test step — run cocotb
    runner.test(
        test_module="Testbench",      # Your Python test file (without the .py extension)
        testcase=target_tests,       # Pass the specific test case(s) here
        hdl_toplevel="ALU",
        hdl_toplevel_lang="verilog",
        verbose=True,
        gui=use_gui,                 # Controls whether Questa opens visually
        test_args=[
            "-l", "transcript.log"   # Saves Questa's terminal output to a log file
        ],
    )

if __name__ == "__main__":
    # Add an argument parser to select tests and modes from the command line
    parser = argparse.ArgumentParser(description="Run ALU Cocotb Tests")
    
    parser.add_argument("--test", type=str, default=None,
                        help="Name of specific test to run (e.g., add_xor_test)")
    
    parser.add_argument("--regression", action="store_true",
                        help="Run the full suite of tests sequentially")
    
    parser.add_argument("--no-gui", action="store_true",
                        help="Run purely in the terminal (disables Questa GUI)")
    
    args = parser.parse_args()

    # Run the function with the parsed arguments (GUI is True by default)
    get_questa_runner(testcase=args.test, run_regression=args.regression, use_gui=not args.no_gui)