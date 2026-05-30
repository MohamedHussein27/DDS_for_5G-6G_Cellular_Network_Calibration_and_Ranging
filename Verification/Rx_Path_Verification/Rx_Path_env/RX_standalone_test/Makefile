# ==============================================================================
# cocotb Makefile: ADI 5G/6G Cellular Network Calibration - RX Path (ISAC)
# ==============================================================================
export COCOTB_RESOLVE_X = ZEROS


# Simulator Defaults
GUI ?= 1
WAVES ?= 1
SIM ?= questa
TOPLEVEL_LANG ?= verilog



# VERILOG SOURCES 
# ─────────────────────────────────────────────────────────
V_FILES = multiplier.v \
          butterfly.v \
          controlunit_2048pts.v \
          controlunit_4096pts.v \
          delayfeedback.v \
          fft__4096_stage.v \
          fft_4096_top.v \
          ifft_2048_stage.v \
          ifft_2048_top.v \
          ifft_4096_stage.v \
          ifft_4096_top.v \
          isac_demux.v \
          bit_reversal.v \
          RAM_RX.v \
          twiddlerom_2048.v \
          twiddlerom_4096pts.v \
          RX_top.v

# $(CURDIR) ensures the simulator finds the files from any directory
VERILOG_SOURCES = $(foreach f, $(V_FILES), $(CURDIR)/$(f))

# ─────────────────────────────────────────────────────────
# Top-level Verilog module name
TOPLEVEL = RX_TOP

# Python test module name (Default fallback, overridden by Makefile.runs)
MODULE = rx_top_test

# Tell the Questa compiler to instrument the RTL for coverage (bcest)
COMPILE_ARGS += +cover=bcest

include $(shell cocotb-config --makefiles)/Makefile.sim

# ─────────────────────────────────────────────────────────
# CLEAN TARGET (Must use a real TAB on the line below!)
clean::
	rm -rf __pycache__ results.xml modelsim.ini transcript vsim.wlf sim_build

