import os
import sys
from pathlib import Path
from cocotb.runner import get_runner


def get_questa_runner():
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
        sources=[alu_file],        # list your design files here
        hdl_toplevel="ALU",        # your top-level module
        always=True,
        clean=True,
        waves=True,                # generate waveform dumps
    )

    # Test step — run cocotb
    runner.test(
        test_module="testbench",
        testcase=None,
        hdl_toplevel="ALU",
        hdl_toplevel_lang="verilog",
        verbose=True,
        gui=True,
        test_args=[
            "-l", "transcript.log"
        ],
    )


if __name__ == "__main__":
    get_questa_runner()
