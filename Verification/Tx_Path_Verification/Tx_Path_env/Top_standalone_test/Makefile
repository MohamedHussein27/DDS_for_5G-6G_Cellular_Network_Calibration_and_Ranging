# ==============================================================================
# PyUVM Makefile for the 4096-point IFFT
# Sponsor  : Analog Devices, Inc. (ADI)
# Project  : DDS for 5G/6G Cellular Network Calibration and Ranging
# ==============================================================================
export COCOTB_RESOLVE_X = ZEROS

# 1. Simulator Configuration (Icarus Verilog)
SIM ?= icarus
TOPLEVEL_LANG ?= verilog

# Enable SystemVerilog parsing for Icarus
COMPILE_ARGS += -g2012

# Enable waveform dumping (generates dump.vcd for GTKWave)
WAVES ?= 1
 
# ------------------------------------------------------------------------------
# RTL source files
# ------------------------------------------------------------------------------
VERILOG_SOURCES += $(PWD)/multiplier.sv
VERILOG_SOURCES += $(PWD)/butterfly.sv
VERILOG_SOURCES += $(PWD)/controlunit_4096.sv
VERILOG_SOURCES += $(PWD)/ofdm_rom.sv
VERILOG_SOURCES += $(PWD)/delayfeedback.sv
VERILOG_SOURCES += $(PWD)/PHASE_ACC.sv
VERILOG_SOURCES += $(PWD)/first_quad_address.sv
VERILOG_SOURCES += $(PWD)/negative_mux.sv
VERILOG_SOURCES += $(PWD)/ifft_4096_stage.sv
VERILOG_SOURCES += $(PWD)/twiddlerom_4096.sv
VERILOG_SOURCES += $(PWD)/ifft_4096_top.sv
VERILOG_SOURCES += $(PWD)/bit_reversal.sv
VERILOG_SOURCES += $(PWD)/config_reg.sv
VERILOG_SOURCES += $(PWD)/DDS.sv
VERILOG_SOURCES += $(PWD)/fft_4096_stage.sv
VERILOG_SOURCES += $(PWD)/fft_4096_top.sv
VERILOG_SOURCES += $(PWD)/isac_mux.sv
VERILOG_SOURCES += $(PWD)/LUT.sv
VERILOG_SOURCES += $(PWD)/TX_TOP.sv

# Top-level hardware module
TOPLEVEL = TX_TOP
 
# ==============================================================================
# Execution Targets
# ==============================================================================
MODULE = stress_test

# 5. UVM Test Selection Hook
# This allows your Python code to catch the specific test name from the terminal.
export UVM_TESTNAME ?= stress_test

# 6. Include standard Cocotb makefiles
include $(shell cocotb-config --makefiles)/Makefile.sim

# ============================================================================
# Cleanup Target
# ============================================================================
.PHONY: clean
clean::
	rm -rf sim_build results.xml __pycache__ transcript.log dump.vcd waves.vcd