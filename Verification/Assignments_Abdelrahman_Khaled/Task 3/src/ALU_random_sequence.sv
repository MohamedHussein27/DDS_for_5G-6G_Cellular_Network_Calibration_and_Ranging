package ALU_random_sequence_pkg;

  import uvm_pkg::*;
  `include "uvm_macros.svh"

  import ALU_seq_item_pkg::*;

  //----------------------------------------------------------------------------
  // ALU Random Sequence
  //----------------------------------------------------------------------------
  // Purpose:
  // - Generates randomized ALU transactions
  // - Randomizes operation (ADD, XOR, AND, OR)
  // - Randomizes operands a and b (4-bit)
  // - Sends transactions to the ALU driver via the sequencer
  //----------------------------------------------------------------------------
  class ALU_random_sequence extends uvm_sequence #(ALU_seq_item);

    // Register sequence with UVM factory
    `uvm_object_utils(ALU_random_sequence)

    // Transaction handle
    ALU_seq_item item;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    function new(string name = "ALU_random_sequence");
      super.new(name);
    endfunction

    //--------------------------------------------------------------------------
    // Body Task
    //--------------------------------------------------------------------------
    // This task runs when the sequence is started on the sequencer
    //--------------------------------------------------------------------------
    virtual task body();

      // Generate a number of random ALU transactions
      repeat (1000) begin

        // Create a new transaction
        item = ALU_seq_item::type_id::create("item");

        // Start sequence item handshake
        start_item(item);

        // Randomize all fields of the transaction
        if (!item.randomize()) begin
          `uvm_error("ALU_RAND_SEQ", "Randomization failed")
        end

        // Finish item and send to driver
        finish_item(item);

      end

    endtask : body

  endclass : ALU_random_sequence

endpackage : ALU_random_sequence_pkg
