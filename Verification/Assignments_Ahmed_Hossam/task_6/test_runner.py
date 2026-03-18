import os
import sys
import argparse
from pathlib import Path
from cocotb.runner import get_runner

def get_questa_runner(seq_name=None, run_regression=False, use_gui=True):
    # Default to Questasim
    SIM = os.getenv("SIM", "questa")
    project_path = Path(__file__).parent.resolve()

    # Make sure Python can find your local modules (environment.py, etc.)
    sys.path.insert(0, str(project_path))

    # Path to your SystemVerilog design files
    # IMPORTANT: Use .sv extensions so vlog recognizes the file type
    sv_sources = [
        project_path.joinpath("DUT.sv.sv").as_posix()
    ]

    # Get the simulator runner
    runner = get_runner(SIM)

    # Build step — compile your DUT
    runner.build(
        verilog_sources=sv_sources,
        hdl_toplevel="ALU",
        always=True,
        build_args=["-sv", "+incdir+" + str(project_path)] # Added -sv flag from reference
    )

    # Determine what to run based on arguments
    if run_regression:
        target_seq = "all"
        print("\n=== RUNNING REGRESSION SUITE (ALL SEQUENCES) ===")
    else:
        target_seq = seq_name if seq_name else "all"
        print(f"\n=== RUNNING SEQUENCE: {target_seq} ===")

    # Test step — run cocotb
    runner.test(
        test_module="testbench",
        hdl_toplevel="ALU",
        hdl_toplevel_lang="verilog",
        verbose=True,
        gui=use_gui,
        waves=True,
        test_args=["-l", "transcript.log"],
        # Injects the variable into os.environ for generator.py
        extra_env={"TEST_SEQ": target_seq} 
    )

if __name__ == "__main__":
    # Add an argument parser similar to your reference code
    parser = argparse.ArgumentParser(description="ALU Cocotb Test Runner")
    
    parser.add_argument("--seq", type=str, default="all",
                        help="Specify the sequence (add, xor, and, or, all)")
    
    parser.add_argument("--regression", action="store_true",
                        help="Run the full suite of sequences")
    
    parser.add_argument("--no-gui", action="store_true",
                        help="Run purely in the terminal (disables Questa GUI)")
    
    args = parser.parse_args()

    # Run the function
    get_questa_runner(
        seq_name=args.seq, 
        run_regression=args.regression, 
        use_gui=not args.no_gui
    )