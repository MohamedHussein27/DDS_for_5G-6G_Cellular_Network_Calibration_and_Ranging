# ==============================================================================
# cocotb Makefile: ADI 5G/6G Cellular Network Calibration - FFT Block
# ==============================================================================
export COCOTB_RESOLVE_X = ZEROS
export PYTHONEXECUTABLE=C:/Users/Moustafa/miniconda3/python.exe

# Simulator Defaults
GUI ?= 1
WAVES ?= 1
SIM ?= questa
TOPLEVEL_LANG ?= verilog

# 1. Point to the FFT SystemVerilog/Verilog RTL files
# 1. Point to the FFT RTL files (Use relative paths for Windows)
# 1. Point to the FFT RTL files (Using direct relative paths for Windows/Questa compatibility)
# VERILOG SOURCES (The "Windows-Proof" Version)
# ─────────────────────────────────────────────────────────

# This forces the path to be a native Windows path that Questa understands
PROJ_DIR := $(shell pwd -W)

# List your files here
# VERILOG SOURCES (The "Universal" Version)
# ─────────────────────────────────────────────────────────

# $(CURDIR) is a built-in Make variable that always holds the absolute path
# to the folder where the Makefile is. No 'pwd' command needed!
V_FILES = twiddlerom_4096pts.v \
          butterfly.v \
          controlunit_4096pts.v \
          delayfeedback.v \
          fft_4096_stage.v \
          multiplier.v \
          fft_4096_top.v

# We use $(CURDIR)/$(f) to make sure the simulator finds the files
# even when running from inside the 'sim_build' folder.
VERILOG_SOURCES = $(foreach f, $(V_FILES), $(CURDIR)/$(f))
# ─────────────────────────────────────────────────────────
# Debug: Uncomment the line below if it fails again to see 
# exactly what paths are being sent to the simulator
# $(info VERILOG_SOURCES is $(VERILOG_SOURCES))
# 2. Top-level Verilog module name
TOPLEVEL = fft_4096_top

# 3. Python test module name (matches test_fft.py)
MODULE = fft_test

# Include standard cocotb makefiles
include $(shell cocotb-config --makefiles)/Makefile.sim

# Clean target
.PHONY: clean
clean::
	rm -rf __pycache__ results.xml modelsim.ini transcript vsim.wlf sim_build