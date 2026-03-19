# Makefile
# ==============================================================================
# Cocotb & Simulator Configuration
# ==============================================================================
export COCOTB_RESOLVE_X = ZEROS

# Defaults
GUI ?= 1
SIM ?= questa
TOPLEVEL_LANG ?= verilog

# Point to the Verilog file in the current directory
VERILOG_SOURCES += $(PWD)/alu.sv

# Match the module name inside alu.sv
TOPLEVEL = ALU

# Match your Python entry-point file name (e.g., test_top.py)
MODULE = test_top

# Include standard Cocotb makefiles
include $(shell cocotb-config --makefiles)/Makefile.sim

# ==============================================================================
# Cleanup Target
# ==============================================================================
.PHONY: clean
clean::
	rm -rf __pycache__ results.xml modelsim.ini transcript vsim.wlf sim_build