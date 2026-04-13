import os
import sys
from pathlib import Path
from cocotb.runner import get_runner


def run():

    SIM = os.getenv("SIM", "questa")
    project_path = Path(__file__).parent.resolve()

    # allow python to find your testbench
    sys.path.insert(0, str(project_path / "src"))

    # DUT file
    alu_file = project_path / "src" / "dut.sv"

    runner = get_runner(SIM)

    # -------------------------------------------------
    # BUILD
    # -------------------------------------------------
    runner.build(
        sources=[alu_file],
        hdl_toplevel="ALU",
        always=True,
        clean=True,
        waves=True,
    )

    # -------------------------------------------------
    # TEST SELECTION (LIKE MAKE TEST=add)
    # -------------------------------------------------
    test = os.getenv("TEST", "regression")

    print(f"\n==============================")
    print(f" RUNNING TEST: {test}")
    print(f"==============================\n")

    # -------------------------------------------------
    # RUN
    # -------------------------------------------------
    runner.test(
        test_module="ALU_test",   # your python file name
        hdl_toplevel="ALU",
        hdl_toplevel_lang="verilog",
        verbose=True,
        gui=True,
        test_args=["-l", "transcript.log"],
        extra_env={
            "UVM_TESTNAME": test
        }
    )


if __name__ == "__main__":
    run()
