import os
import argparse
import subprocess  # <--- Imported to run terminal commands from Python
from pathlib import Path
from cocotb.runner import get_runner

def test_alu_runner():
    # 1. Set up the argument parser
    parser = argparse.ArgumentParser(description="Run ALU UVM simulation.")
    parser.add_argument(
        "-t", "--test", 
        type=str, 
        default="test_all",
        help="Name of the UVM test class to run (e.g., test_addition, test_xor)"
    )
    # --- NEW: Add a flag for waveforms ---
    parser.add_argument(
        "-w", "--wave", 
        action="store_true", # Sets args.wave to True if you type -w
        help="Open Questa waveform viewer after the simulation finishes"
    )
    args = parser.parse_args()

    # 2. Setup the runner and paths
    sim = os.getenv("SIM", "questa")
    project_path = Path(__file__).resolve().parent
    sv_sources = [project_path.joinpath("DUT.sv.sv").as_posix()]

    runner = get_runner(sim)
    
    # 3. Build the design
    runner.build(
        verilog_sources=sv_sources,
        hdl_toplevel="ALU",
        always=True,
        build_args=["-sv"] 
    )
    
    # 4. Run the simulation
    runner.test(
        hdl_toplevel="ALU",
        test_module="test", 
        extra_env={"UVM_TESTNAME": args.test},
        waves=True 
    )

    # 5. --- NEW: Launch Waveforms if requested ---
    if args.wave:
        print("\n>>> Launching Questa Waveform Viewer... <<<")
        wlf_path = project_path / "sim_build" / "vsim.wlf"
        
        # Check if the file actually generated successfully
        if wlf_path.exists():
            subprocess.run(["vsim", "-view", str(wlf_path)])
        else:
            print("Error: Waveform file not found. The simulation might have crashed.")

if __name__ == "__main__":
    test_alu_runner()