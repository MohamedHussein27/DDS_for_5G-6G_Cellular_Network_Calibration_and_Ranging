package ALU_sequencer_pkg;

  // Import UVM base package and macros
  import uvm_pkg::*;
  `include "uvm_macros.svh"

  // Import ALU sequence item definition
  import ALU_seq_item_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Sequencer
  //----------------------------------------------------------------------------
  // Responsibilities:
  // 1) Acts as an arbitrator between sequences and the driver
  // 2) Provides sequence items to the driver upon request
  // 3) Controls the flow of stimulus
  //
  // Note:
  // - The sequencer itself does NOT generate transactions
  // - Sequences run on the sequencer and produce ALU_seq_item objects
  //----------------------------------------------------------------------------
  class ALU_sequencer extends uvm_sequencer #(ALU_seq_item);

    // Register sequencer with UVM factory
    `uvm_component_utils(ALU_sequencer)

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_sequencer", uvm_component parent = null);
      super.new(name, parent);
    endfunction

  endclass : ALU_sequencer

endpackage : ALU_sequencer_pkg
