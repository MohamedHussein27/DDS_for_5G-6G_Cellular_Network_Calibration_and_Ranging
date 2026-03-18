# ==============================================================================
# Multi-File PyUVM Makefile for Icarus Verilog
# ==============================================================================

SIM ?= icarus
TOPLEVEL_LANG ?= verilog

# Enable SystemVerilog syntax in Icarus
COMPILE_ARGS += -g2012

# Define the source files (Hardware RTL only)
VERILOG_SOURCES += $(PWD)/DUT.sv

# Define the Top-Level hardware module
TOPLEVEL = ALU

# ==============================================================================
# Execution Targets
# ==============================================================================

# Default target: The Regression Suite
# By separating the file names with commas (NO SPACES), Cocotb will run them all sequentially!
regression:
	$(MAKE) clean_sim
	$(MAKE) add_xor
	$(MAKE) clean_sim
	$(MAKE) and_or
	$(MAKE) clean_sim
	$(MAKE) rand

# Individual Test Targets
add_xor:
	$(MAKE) sim MODULE=add_xor_test

and_or:
	$(MAKE) sim MODULE=and_or_test

rand:
	$(MAKE) sim MODULE=rand_test

# Clean up simulation artifacts
clean_sim:
	$(MAKE) clean
	rm -rf sim_build results.xml __pycache__ transcript.log dump.vcd

# ==============================================================================
# Include standard Cocotb machinery
# ==============================================================================
include $(shell cocotb-config --makefiles)/Makefile.sim