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
# 2. RTL Source Files
# ------------------------------------------------------------------------------
VERILOG_SOURCES += $(PWD)/bit_reversal.v
VERILOG_SOURCES += $(PWD)/butterfly.v
VERILOG_SOURCES += $(PWD)/config_reg.v
VERILOG_SOURCES += $(PWD)/controlunit_2048pts.v
VERILOG_SOURCES += $(PWD)/controlunit_4096pts.v
VERILOG_SOURCES += $(PWD)/DDS.v
VERILOG_SOURCES += $(PWD)/delayfeedback.v
VERILOG_SOURCES += $(PWD)/fft__4096_stage.v
VERILOG_SOURCES += $(PWD)/fft_4096_stage.v
VERILOG_SOURCES += $(PWD)/fft_4096_top.v
VERILOG_SOURCES += $(PWD)/first_quad_address.v
VERILOG_SOURCES += $(PWD)/ifft_2048_stage.v
VERILOG_SOURCES += $(PWD)/ifft_2048_top.v
VERILOG_SOURCES += $(PWD)/ifft_4096_stage.v
VERILOG_SOURCES += $(PWD)/ifft_4096_top.v
VERILOG_SOURCES += $(PWD)/isac_demux.v
VERILOG_SOURCES += $(PWD)/isac_mux.v
VERILOG_SOURCES += $(PWD)/LUT.v
VERILOG_SOURCES += $(PWD)/multiplier.v
VERILOG_SOURCES += $(PWD)/multiplier_RX.v
VERILOG_SOURCES += $(PWD)/negative_mux.v
VERILOG_SOURCES += $(PWD)/ofdm_rom.v
VERILOG_SOURCES += $(PWD)/PHASE_ACC.v
VERILOG_SOURCES += $(PWD)/RAM_RX.v
VERILOG_SOURCES += $(PWD)/RX_top.v
VERILOG_SOURCES += $(PWD)/SYSTEM_TOP.v
VERILOG_SOURCES += $(PWD)/tb_fft_4096.v
VERILOG_SOURCES += $(PWD)/twiddlerom_2048.v
VERILOG_SOURCES += $(PWD)/twiddlerom_4096pts.v
VERILOG_SOURCES += $(PWD)/TX_TOP.v

# 3. Top-level Verilog module name
TOPLEVEL = SYSTEM_TOP

# ==============================================================================
# Default Execution Targets (Fallbacks)
# ==============================================================================
# If you just run 'make sim', it will default to the normal_test
MODULE ?= system_top_normal_test
export UVM_TESTNAME ?= system_top_normal_test

# 6. Include standard Cocotb makefiles
include $(shell cocotb-config --makefiles)/Makefile.sim

# ==============================================================================
# Custom Test Functions (Easy Terminal Commands)
# ==============================================================================
.PHONY: normal clean

# Run the Normal Test
normal:
	$(MAKE) sim MODULE=system_top_normal_test UVM_TESTNAME=system_top_normal_test

# ============================================================================
# Cleanup Target
# ============================================================================
clean::
	rm -rf sim_build results.xml __pycache__ transcript.log dump.vcd waves.vcd
