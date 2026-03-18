# ==============================================================================
# PyUVM Makefile for Icarus Verilog (Decorator Mode)
# ==============================================================================

SIM ?= icarus
TOPLEVEL_LANG ?= verilog

# Enable SystemVerilog syntax in Icarus
COMPILE_ARGS += -g2012

# Define the source files (Hardware RTL only)
VERILOG_SOURCES += $(PWD)/DUT.sv

# Define the Top-Level hardware module
TOPLEVEL = ALU

# Define the Python test file
MODULE = test

# ==============================================================================
# Execution Targets
# ==============================================================================

# Default target runs all @pyuvm.test() decorated classes in test.py
all: sim

clean_sim:
	$(MAKE) clean
	rm -rf sim_build results.xml __pycache__ transcript.log

# ==============================================================================
# Include standard Cocotb machinery
# ==============================================================================
include $(shell cocotb-config --makefiles)/Makefile.sim