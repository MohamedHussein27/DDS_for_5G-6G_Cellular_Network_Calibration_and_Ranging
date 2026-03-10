# Makefile
export COCOTB_RESOLVE_X = ZEROS

# defaults
GUI=1
WAVES=1
SIM ?= questa
TOPLEVEL_LANG ?= verilog

# 1.  Point to the new Verilog file in the current directory
VERILOG_SOURCES += $(PWD)/alu.sv

# 2.  Change TOPLEVEL to match the module name inside alu.sv
TOPLEVEL = ALU

# MODULE is already correct (matches test.py)
MODULE = test

include $(shell cocotb-config --makefiles)/Makefile.sim

# Clean target
.PHONY: clean
clean::
	rm -rf __pycache__ results.xml modelsim.ini transcript vsim.wlf

