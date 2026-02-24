# Makefile
export COCOTB_RESOLVE_X = ZEROS

# defaults
GUI=1
WAVES=1
SIM ?= questa
TOPLEVEL_LANG ?= verilog

# 1.  Point to the new Verilog file in the current directory
VERILOG_SOURCES += $(PWD)/simple_dut.sv

# 2.  Change TOPLEVEL to match the module name inside simple_dut.sv
TOPLEVEL = simple_dut

# MODULE is already correct (matches test_top.py)
MODULE = test_top

include $(shell cocotb-config --makefiles)/Makefile.sim

# Clean target
.PHONY: clean
clean::
	rm -rf __pycache__ results.xml modelsim.ini transcript vsim.wlf

