# ============================================================================
# Full Verification Makefile (pyuvm + cocotb + Questa)
# ============================================================================

# Handle undefined signals (from your old makefile)
export COCOTB_RESOLVE_X = ZEROS

# 1. Simulator Configuration (Using your working setup)
GUI ?= 1
SIM ?= questa
TOPLEVEL_LANG ?= verilog

# 2. Design Under Test (DUT) Sources
# Point to your actual SystemVerilog file (from your screenshot)
# 2. Design Under Test (DUT) Sources
VERILOG_SOURCES += DDS.sv \
                   PHASE_ACC.sv \
                   first_quad_address.sv \
                   LUT.sv \
                   negative_mux.sv \
                   bit_reversal.sv \
                   butterfly.sv \
                   config_reg.sv \
                   controlunit_4096.sv \
                   delayfeedback.sv \
                   fft_4096_stage.sv \
                   fft_4096_top.sv \
                   ifft_4096_stage.sv \
                     ifft_4096_top.sv \
                     isac_mux.sv \
                     multiplier.sv \
                     ofdm_rom.sv \
                     twiddlerom_4096.sv \
                     TX_TOP.sv 

# 3. Top-Level Module Name
# This MUST match the exact module name inside DUT.sv (Case-Sensitive!)
TOPLEVEL = TX_TOP

# 4. Python Test Module
# The name of your Python file containing the @cocotb.test() function (WITHOUT the .py)
# NOTE: Replace "test_alu" with whatever your python file is actually named!
MODULE = reset_test

# 5. UVM Test Selection Hook
# This allows your Python code to catch the specific test name from the terminal.
export UVM_TESTNAME ?= reset_test_only

# 6. Include standard Cocotb makefiles
include $(shell cocotb-config --makefiles)/Makefile.sim

# ============================================================================
# Cleanup Target
# ============================================================================
.PHONY: clean
clean::
	rm -rf __pycache__ results.xml transcript sim_build dds_coverage.xml modelsim.ini vsim.wlf