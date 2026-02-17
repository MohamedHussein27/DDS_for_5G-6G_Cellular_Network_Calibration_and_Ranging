package ALU_seq_item_pkg;

  import uvm_pkg::*;
  `include "uvm_macros.svh"

  //----------------------------------------------------------------------------
  // ALU Sequence Item
  //----------------------------------------------------------------------------
  // Represents one ALU transaction:
  //   - Inputs  : a, b, op
  //   - Outputs : out, c
  //----------------------------------------------------------------------------
  class ALU_seq_item extends uvm_sequence_item;

    //--------------------------------------------------------------------------
    // UVM Factory Registration & Field Automation
    //--------------------------------------------------------------------------
    `uvm_object_utils(ALU_seq_item);

    //--------------------------------------------------------------------------
    // ALU Inputs (Randomized)
    //--------------------------------------------------------------------------
    rand logic [3:0] a;          // Operand A
    rand logic [3:0] b;          // Operand B
    rand logic [1:0] op;         // Operation selector

    //--------------------------------------------------------------------------
    // ALU Outputs (Sampled by monitor)
    //--------------------------------------------------------------------------
    logic [3:0] out;             // ALU result
    logic       c;               // Carry-out

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_seq_item");
      super.new(name);
    endfunction

    //--------------------------------------------------------------------------
    // Convert to string (debug-friendly)
    //--------------------------------------------------------------------------
    function string convert2string();
      return $sformatf(
        "ALU_TXN : op=%0d a=0x%0h b=0x%0h | out=0x%0h c=%0b",
        op, a, b, out, c
      );
    endfunction

  endclass : ALU_seq_item

endpackage : ALU_seq_item_pkg
