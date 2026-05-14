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
VERILOG_SOURCES += $(PWD)/multiplier.v
VERILOG_SOURCES += $(PWD)/butterfly.v
VERILOG_SOURCES += $(PWD)/controlunit_4096pts.v
VERILOG_SOURCES += $(PWD)/delayfeedback.v
VERILOG_SOURCES += $(PWD)/ifft_4096_stage.v
VERILOG_SOURCES += $(PWD)/twiddlerom_4096pts.v
VERILOG_SOURCES += $(PWD)/ifft_4096_top.v

# Top-level hardware module
TOPLEVEL = ifft_4096_top
 
# ==============================================================================
# Execution Targets
# ==============================================================================
 
# Default target: full regression
all: regression
 
# ------------------------------------------------------------------------------
# Regression — runs every test in order, cleans between each run
# ------------------------------------------------------------------------------
# Full regression — all sequences in one simulation run
regression:
	$(MAKE) clean_sim
	$(MAKE) reset
	$(MAKE) clean_sim
	$(MAKE) allzero
	$(MAKE) clean_sim
	$(MAKE) impulse
	$(MAKE) clean_sim
	$(MAKE) dc
	$(MAKE) clean_sim
	$(MAKE) singletone
	$(MAKE) clean_sim
	$(MAKE) twotone
	$(MAKE) clean_sim
	$(MAKE) functional
	$(MAKE) clean_sim
	$(MAKE) maxpositive
	$(MAKE) clean_sim
	$(MAKE) maxnegative
	$(MAKE) clean_sim
	$(MAKE) halfnyquist
	$(MAKE) clean_sim
	$(MAKE) alternating
	$(MAKE) clean_sim
	$(MAKE) burst
	$(MAKE) clean_sim
	$(MAKE) ramp
	$(MAKE) clean_sim
	$(MAKE) backtoback
	$(MAKE) clean_sim
	$(MAKE) random
	$(MAKE) clean_sim
 
# ------------------------------------------------------------------------------
# Individual test targets
# ------------------------------------------------------------------------------
 
# Reset only — rst_n flushes registers, valid_out stays 0
reset:
	$(MAKE) sim MODULE=ifft_reset_test

# allzero test
allzero:
	$(MAKE) sim MODULE=ifft_allzero_test

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

maxpositive:
	$(MAKE) sim MODULE=ifft_maxpositive_test

maxnegative:
	$(MAKE) sim MODULE=ifft_maxnegative_test

halfnyquist:
	$(MAKE) sim MODULE=ifft_halfnyquist_test

alternating:
	$(MAKE) sim MODULE=ifft_alternating_test

burst:
	$(MAKE) sim MODULE=ifft_burst_test

ramp:
	$(MAKE) sim MODULE=ifft_ramp_test

backtoback:
	$(MAKE) sim MODULE=ifft_backtoback_test

random:
	$(MAKE) sim MODULE=ifft_random_test



# Frequency-axis sweep — twiddle ROM at bins 1, 100, 2047 + Nyquist
frequency:
	$(MAKE) sim MODULE=ifft_frequency_test
 
# Datapath limits — saturation, carry/borrow, full dynamic range
stress:
	$(MAKE) sim MODULE=ifft_stress_test
 
 
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