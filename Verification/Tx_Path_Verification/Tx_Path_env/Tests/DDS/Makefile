# ==============================================================================
# cocotb Makefile: SPI Master Verification (pyuvm + cocotb + Questa)
# ==============================================================================

# Handle undefined signals as zeros (prevents X-propagation surprises)
export COCOTB_RESOLVE_X = ZEROS



# ==============================================================================
# 1. Simulator Configuration
# ==============================================================================
GUI   ?= 1
WAVES ?= 1
SIM   ?= questa
TOPLEVEL_LANG ?= verilog

# ==============================================================================
# 2. DUT Source Files
#    PROJ_DIR forces a native Windows path that Questa understands when
#    running inside the Git-Bash / MSYS2 shell.
#    $(CURDIR)/$(f) ensures the simulator finds files even when it spawns
#    from inside the sim_build/ subdirectory.
# ==============================================================================
PROJ_DIR := $(shell pwd -W)

V_FILES = DDS.v \
          PHASE_ACC.v \
          first_quad_address.v \
          LUT.v \
          negative_mux.v

# Expand every filename to its absolute path
VERILOG_SOURCES = $(foreach f, $(V_FILES), $(CURDIR)/$(f))

# ==============================================================================
# 3. Top-Level Module
# ==============================================================================
TOPLEVEL = dds_top

# ==============================================================================
# 4. Coverage Instrumentation
#    +cover=bcest tells Questa to instrument for:
#      b = branch  c = condition  e = expression  s = statement  t = toggle
# ==============================================================================
COMPILE_ARGS += +cover=bcest

# ==============================================================================
# Include the standard cocotb simulation rules
# (defines the 'sim' target, compiles the DUT, launches the simulator)
# ==============================================================================
include $(shell cocotb-config --makefiles)/Makefile.sim

# ==============================================================================
# 5. Execution Targets  (individual runs live in Makefile.runs)
# ==============================================================================
include Makefile.runs

# ==============================================================================
# 6. Regression — runs every test module in a single simulator invocation.
#    Comma-separated MODULE list lets cocotb discover and run all tests
#    without a separate clean between each one.
# ==============================================================================
regression:
	$(MAKE) sim \
	  MODULE=dds_reset_test,dds_random_test,dds_chirpsweep_test,\
dds_cyclestress_test,dds_single_tone_test,dds_fmcw_test,\
dds_cov_test,dds_reset_recovery_test,dds_psk_modulation_test,\
dds_fft_boundry_test,dds_dc_hold_test \
	  GUI=1

# ==============================================================================
# 7. Cleanup Targets
# ==============================================================================
.PHONY: clean clean_sim all regression

clean_sim:
	rm -rf __pycache__ results.xml transcript sim_build \
	       dds_coverage.xml modelsim.ini vsim.wlf

clean:: clean_sim
