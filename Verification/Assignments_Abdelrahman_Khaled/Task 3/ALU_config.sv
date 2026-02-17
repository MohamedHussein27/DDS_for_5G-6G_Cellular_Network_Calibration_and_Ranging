package ALU_config_pkg;

  // Import UVM base package
  import uvm_pkg::*;
  `include "uvm_macros.svh"

  //----------------------------------------------------------------------------
  // ALU Configuration Object
  //----------------------------------------------------------------------------
  // Purpose:
  // - Hold configuration information for ALU agent components
  // - Passed through uvm_config_db
  // - Shared between driver, monitor, and sequencer if needed
  //----------------------------------------------------------------------------
  class ALU_config extends uvm_object;

    // Register config object with UVM factory
    `uvm_object_utils(ALU_config)

    // Virtual interface to connect TB components to DUT
    virtual ALU_if vif;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_config");
      super.new(name);
    endfunction

  endclass : ALU_config

endpackage : ALU_config_pkg
