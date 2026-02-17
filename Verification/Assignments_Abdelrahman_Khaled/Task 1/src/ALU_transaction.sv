// Package that contains the ALU transaction class
package ALU_transaction_pkg;

  // Transaction class representing one ALU operation
  class ALU_transaction;

    // Randomized input operands
    randc logic [3:0] a, b;      // ALU input operands

    // Randomized operation selector
    randc logic [1:0] op;        // ALU operation code

    // DUT outputs (not randomized)
    logic [3:0] out;            // ALU result
    logic       c;              // Carry-out flag


  endclass : ALU_transaction

endpackage : ALU_transaction_pkg
