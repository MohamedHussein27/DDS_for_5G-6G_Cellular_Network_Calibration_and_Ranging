# Define the simulator (Change this to whatever you use: icarus, questa, verilator, vcs)
SIM := icarus

# Define the source files (Hardware RTL only)
VERILOG_SOURCES += $(PWD)/DUT.sv

# Define the Top-Level hardware module
TOPLEVEL = ALU

# Define the Python test module (the name of your Python file without the .py)
MODULE = test

# ------------------------------------------------------------------------------
# PyUVM Magic: Pass the class name to the simulator
# NOTE: If you use Questa/VCS, change PLUSARGS to EXTRA_ARGS
# ------------------------------------------------------------------------------
TEST_NAME ?= test_all
PLUSARGS += +UVM_TESTNAME=$(TEST_NAME)

# ==============================================================================
# Execution Targets
# ==============================================================================

# Default target: The Regression Suite
regression:
	$(MAKE) clean_sim
	$(MAKE) add
	$(MAKE) clean_sim
	$(MAKE) xor
	$(MAKE) clean_sim
	$(MAKE) and
	$(MAKE) clean_sim
	$(MAKE) or
	$(MAKE) clean_sim
	$(MAKE) all

# Individual Test Targets (These overwrite the TEST_NAME variable dynamically)
add:
	$(MAKE) sim TEST_NAME=test_addition

xor:
	$(MAKE) sim TEST_NAME=test_xor

and:
	$(MAKE) sim TEST_NAME=test_and

or:
	$(MAKE) sim TEST_NAME=test_or

all:
	$(MAKE) sim TEST_NAME=test_all

# Clean up simulation artifacts (Added alu_coverage.xml to the cleanup list!)
clean_sim:
	$(MAKE) clean
	rm -rf sim_build results.xml __pycache__ transcript.log dump.vcd alu_coverage.xml

# ==============================================================================
# Include standard Cocotb machinery
# ==============================================================================
include $(shell cocotb-config --makefiles)/Makefile.sim