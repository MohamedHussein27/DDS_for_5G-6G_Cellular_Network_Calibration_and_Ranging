package ALU_driver_pkg;

  // Import UVM base library
  import uvm_pkg::*;
  `include "uvm_macros.svh"

  // Import transaction and configuration packages
  import ALU_seq_item_pkg::*;
  import ALU_config_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Driver
  //----------------------------------------------------------------------------
  // Responsibilities:
  // 1) Receive sequence items from the sequencer
  // 2) Drive DUT input signals through virtual interface
  // 3) Synchronize with clock
  // 4) Notify sequencer when transaction is completed
  //----------------------------------------------------------------------------
  class ALU_driver extends uvm_driver #(ALU_seq_item);

    // Register this component with the UVM factory
    `uvm_component_utils(ALU_driver)

    // Handle to configuration object
    ALU_config cfg;

    // Virtual interface (obtained from config object)
    virtual ALU_if vif;

    // Transaction handle used for driving
    ALU_seq_item item_dr;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    // Initializes the driver component
    //--------------------------------------------------------------------------
    function new(string name = "ALU_driver", uvm_component parent = null);
      super.new(name, parent);
    endfunction

    //--------------------------------------------------------------------------
    // Run Phase
    //--------------------------------------------------------------------------
    // Main execution phase of the driver.
    // Repeats forever:
    //   1) Get next transaction from sequencer
    //   2) Drive signals to DUT
    //   3) Wait for synchronization (clock edge)
    //   4) Inform sequencer that transaction is done
    //--------------------------------------------------------------------------
    task run_phase(uvm_phase phase);
      super.run_phase(phase);

      forever begin

        // Block until sequencer provides next transaction
        item_dr = ALU_seq_item::type_id::create("item_dr");
        seq_item_port.get_next_item(item_dr);

        // Drive ALU input signals onto DUT interface
        vif.a  = item_dr.a;
        vif.b  = item_dr.b;
        vif.op = item_dr.op;

        // Synchronize driving with clock (negedge example)
        @(negedge vif.clk);

        // Notify sequencer that this transaction is completed
        seq_item_port.item_done();
        `uvm_info("run_phase", item_dr.convert2string(), UVM_HIGH);
      end

    endtask

  endclass : ALU_driver

endpackage : ALU_driver_pkg
