package ALU_agent_pkg;

  // Import UVM base package and macros
  import uvm_pkg::*;
  `include "uvm_macros.svh"

  // Import ALU components
  import ALU_seq_item_pkg::*;
  import ALU_sequencer_pkg::*;
  import ALU_driver_pkg::*;
  import ALU_monitor_pkg::*;
  import ALU_config_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Agent
  //----------------------------------------------------------------------------
  // Responsibilities:
  // 1) Create and connect sequencer, driver, and monitor
  // 2) Distribute configuration to sub-components
  // 3) Provide a single analysis port to the environment
  //
  // The agent can be configured as:
  // - Active  : sequencer + driver + monitor
  // - Passive: monitor only
  //----------------------------------------------------------------------------
  class ALU_agent extends uvm_agent;

    // Register agent with UVM factory
    `uvm_component_utils(ALU_agent)

    // Agent sub-components
    ALU_sequencer sr_ag;   // Sequencer
    ALU_driver    dr_ag;   // Driver
    ALU_monitor   mon_ag;  // Monitor

    // Configuration object
    ALU_config config_ag;

    // Analysis port to forward monitor transactions
    uvm_analysis_port #(ALU_seq_item) agent_ap;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_agent", uvm_component parent = null);
      super.new(name, parent);
    endfunction


    //--------------------------------------------------------------------------
    // Build Phase
    //--------------------------------------------------------------------------
    // - Get configuration object from config_db
    // - Create agent components
    //--------------------------------------------------------------------------
    function void build_phase(uvm_phase phase);
      super.build_phase(phase);

      // Retrieve configuration object
      if (!uvm_config_db#(ALU_config)::get(this, "", "CFG", config_ag))
        `uvm_fatal("AGENT", "Unable to get ALU_config object")

      // Create monitor (always created)
      mon_ag = ALU_monitor::type_id::create("mon_ag", this);

      // Create sequencer and driver only if agent is active
      sr_ag = ALU_sequencer::type_id::create("sr_ag", this);
      dr_ag = ALU_driver   ::type_id::create("dr_ag", this);


      // Create agent-level analysis port
      agent_ap = new("agent_ap", this);
    endfunction


    //--------------------------------------------------------------------------
    // Connect Phase
    //--------------------------------------------------------------------------
    // - Connect driver to sequencer
    // - Connect virtual interface
    // - Forward monitor analysis port
    //--------------------------------------------------------------------------
    function void connect_phase(uvm_phase phase);
      super.connect_phase(phase);

      // Connect virtual interface to driver and monitor
      mon_ag.vif = config_ag.vif;
      dr_ag.vif = config_ag.vif;

    // Connect sequencer to driver
      dr_ag.seq_item_port.connect(sr_ag.seq_item_export);

      // Forward monitor transactions to agent analysis port
      mon_ag.monitor_ap.connect(agent_ap);
    endfunction

  endclass : ALU_agent

endpackage : ALU_agent_pkg
