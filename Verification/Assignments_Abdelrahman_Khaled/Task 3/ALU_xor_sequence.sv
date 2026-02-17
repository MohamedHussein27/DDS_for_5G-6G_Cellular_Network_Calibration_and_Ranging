package ALU_xor_sequence_pkg;

  import uvm_pkg::*;
  `include "uvm_macros.svh"

  import ALU_seq_item_pkg::*;

  //----------------------------------------------------------------------------
  // ALU XOR Sequence
  //----------------------------------------------------------------------------
  // Purpose:
  // - Generates only XOR operations
  // - Iterates over all possible values of a and b (4-bit)
  // - Sends transactions to the ALU driver via the sequencer
  //----------------------------------------------------------------------------
  class ALU_xor_sequence extends uvm_sequence #(ALU_seq_item);

    // Register sequence with UVM factory
    `uvm_object_utils(ALU_xor_sequence)

    // Transaction handle
    ALU_seq_item item;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_xor_sequence");
      super.new(name);
    endfunction

    //--------------------------------------------------------------------------
    // Body Task
    //--------------------------------------------------------------------------
    // This task is executed when the sequence is started on a sequencer
    //--------------------------------------------------------------------------
    virtual task body();

      // Loop over all possible values of operand a (4-bit)
      for (int a = 0; a < 16; a++) begin

        // Loop over all possible values of operand b (4-bit)
        for (int b = 0; b < 16; b++) begin

          // Create a new transaction
          item = ALU_seq_item::type_id::create("item");

          // Start the item (handshake with sequencer/driver)
          start_item(item);

          // Assign transaction fields
          item.op = 2'b01;        // xor operation
          item.a  = a[3:0];
          item.b  = b[3:0];

          // Finish the item (send to driver)
          finish_item(item);

        end
      end

    endtask : body

  endclass : ALU_xor_sequence

endpackage : ALU_xor_sequence_pkg
