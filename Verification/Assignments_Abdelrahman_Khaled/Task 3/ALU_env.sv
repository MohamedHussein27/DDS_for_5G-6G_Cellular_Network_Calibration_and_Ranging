package ALU_env_pkg;

  import uvm_pkg::*;
  `include "uvm_macros.svh"

  import ALU_agent_pkg::*;
  import ALU_scoreboard_pkg::*;
  import ALU_coverage_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Environment
  //----------------------------------------------------------------------------
  // Top-level environment container
  // Instantiates:
  //   - Agent  (driver + sequencer + monitor)
  //   - Scoreboard
  //   - Coverage
  // Connects monitor analysis port to scoreboard & coverage
  //----------------------------------------------------------------------------
  class ALU_env extends uvm_env;

    `uvm_component_utils(ALU_env)

    // Components inside environment
    ALU_agent      ag_env;
    ALU_scoreboard sb_env;
    ALU_coverage   cov_env;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_env", uvm_component parent = null);
      super.new(name, parent);
    endfunction

    //--------------------------------------------------------------------------
    // Build Phase
    // Create all sub-components
    //--------------------------------------------------------------------------
    function void build_phase(uvm_phase phase);
      super.build_phase(phase);

      ag_env  = ALU_agent::type_id::create("ag_env",  this);
      sb_env  = ALU_scoreboard::type_id::create("sb_env",  this);
      cov_env = ALU_coverage::type_id::create("cov_env", this);

    endfunction

    //--------------------------------------------------------------------------
    // Connect Phase
    // Connect analysis ports between components
    //--------------------------------------------------------------------------
    function void connect_phase(uvm_phase phase);
      super.connect_phase(phase);

      // Monitor (inside agent) sends transactions through agent_ap
      // Connect them to:
      //   1) Scoreboard
      //   2) Coverage
      ag_env.agent_ap.connect(sb_env.sb_export);
      ag_env.agent_ap.connect(cov_env.cov_export);

    endfunction

  endclass

endpackage
