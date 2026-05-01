# ==============================================================================
# PyUVM Makefile for the 4096-point IFFT
# Sponsor  : Analog Devices, Inc. (ADI)
# Project  : DDS for 5G/6G Cellular Network Calibration and Ranging
# ==============================================================================
 
SIM            ?= icarus
TOPLEVEL_LANG  ?= verilog
 
# Enable SystemVerilog / Verilog-2012 syntax in Icarus
COMPILE_ARGS   += -g2012
 
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
 
# Default target: full regression
all: regression
 
# ------------------------------------------------------------------------------
# Regression — runs every test in order, cleans between each run
# ------------------------------------------------------------------------------
regression:
	$(MAKE) clean_sim
	$(MAKE) reset
	$(MAKE) clean_sim
	$(MAKE) functional
	$(MAKE) clean_sim
	$(MAKE) frequency
	$(MAKE) clean_sim
	$(MAKE) stress
	$(MAKE) clean_sim
	$(MAKE) random
	$(MAKE) clean_sim
	$(MAKE) pipeline
	$(MAKE) clean_sim
	$(MAKE) reset_recovery
	$(MAKE) clean_sim
	$(MAKE) full_regression
 
# ------------------------------------------------------------------------------
# Individual test targets
# ------------------------------------------------------------------------------
 
# Reset only — rst_n flushes registers, valid_out stays 0
reset:
	$(MAKE) sim MODULE=reset_test

# allzero test
stress:
	$(MAKE) sim MODULE=stress_test

# impulse test
impulse:
	$(MAKE) sim MODULE=ifft_impulse_test

# DC test
dc:
	$(MAKE) sim MODULE=ifft_dc_test

# single-tone test
singletone:
	$(MAKE) sim MODULE=ifft_singletone_test

# two-tone test
twotone:
	$(MAKE) sim MODULE=ifft_twotone_test
 
# Core IFFT math — zero, impulse, DC, single-tone, two-tone linearity
functional:
	$(MAKE) sim MODULE=ifft_functional_test
 



# Frequency-axis sweep — twiddle ROM at bins 1, 100, 2047 + Nyquist
frequency:
	$(MAKE) sim MODULE=ifft_frequency_test
 
# Datapath limits — saturation, carry/borrow, full dynamic range
ifft_stress:
	$(MAKE) sim MODULE=ifft_stress_test
 
# Broad coverage — two random frames with fixed seeds
random:
	$(MAKE) sim MODULE=ifft_random_test
 
# Pipeline timing — back-to-back frames and partial burst
pipeline:
	$(MAKE) sim MODULE=ifft_pipeline_test
 
# Reset recovery — mid-frame reset then clean output verification
reset_recovery:
	$(MAKE) sim MODULE=ifft_reset_recovery_test
 
# Full regression — all 16 sequences in one simulation run
full_regression:
	$(MAKE) sim MODULE=ifft_regression_test
 
# ------------------------------------------------------------------------------
# Clean up simulation artifacts
# ------------------------------------------------------------------------------
clean_sim:
	$(MAKE) clean
	rm -rf sim_build results.xml __pycache__ transcript.log dump.vcd
 
# ==============================================================================
# Include standard Cocotb machinery
# ==============================================================================
include $(shell cocotb-config --makefiles)/Makefile.sim