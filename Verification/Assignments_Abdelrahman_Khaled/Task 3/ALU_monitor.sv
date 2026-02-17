`timescale 1ns/10ps

package ALU_monitor_pkg;

  // Import UVM base package and macros
  import uvm_pkg::*;
  `include "uvm_macros.svh"

  // Import ALU sequence item definition
  import ALU_seq_item_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Monitor
  //----------------------------------------------------------------------------
  // Responsibilities:
  // 1) Passively sample DUT interface signals
  // 2) Package sampled values into ALU_seq_item transactions
  // 3) Send transactions to subscribers (scoreboard, coverage, etc.)
  //----------------------------------------------------------------------------
  class ALU_monitor extends uvm_monitor;

    // Register component with UVM factory
    `uvm_component_utils(ALU_monitor)

    // Transaction used to store sampled data
    ALU_seq_item item_monitor;

    // Analysis port to broadcast sampled transactions
    uvm_analysis_port #(ALU_seq_item) monitor_ap;

    // Virtual interface handle to DUT
    virtual ALU_if vif;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_monitor", uvm_component parent = null);
      super.new(name, parent);
    endfunction


    //--------------------------------------------------------------------------
    // Build Phase
    //--------------------------------------------------------------------------
    // Create analysis port
    //--------------------------------------------------------------------------
    function void build_phase(uvm_phase phase);
      super.build_phase(phase);

      // Initialize analysis port
      monitor_ap = new("monitor_ap", this);
    endfunction


    //--------------------------------------------------------------------------
    // Run Phase
    //--------------------------------------------------------------------------
    // Continuously samples ALU signals on the clock edge and
    // sends the sampled transaction through the analysis port
    //--------------------------------------------------------------------------
    task run_phase(uvm_phase phase);
      super.run_phase(phase);

      forever begin

        // Create a new transaction object for each sample
        item_monitor = ALU_seq_item::type_id::create("item_monitor");

        // Synchronize sampling with DUT clock
        @(posedge vif.clk);
        #0.01; // Small delay to avoid race with driver

        // Sample ALU inputs and outputs
        item_monitor.a   = vif.a;
        item_monitor.b   = vif.b;
        item_monitor.op  = vif.op;
        item_monitor.out = vif.out;
        item_monitor.c   = vif.c;

        // Send sampled transaction to subscribers
        monitor_ap.write(item_monitor);

        // Optional debug print
        `uvm_info("ALU_MON", item_monitor.convert2string(), UVM_HIGH)

      end
    endtask

  endclass : ALU_monitor

endpackage : ALU_monitor_pkg
