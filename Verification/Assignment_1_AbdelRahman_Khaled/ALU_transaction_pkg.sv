// Package that contains the ALU transaction class
package ALU_transaction_pkg;

  // Transaction class representing one ALU operation
  class ALU_transaction;

    // Randomized input operands
    rand logic [3:0] a, b;      // ALU input operands

    // Randomized operation selector
    rand logic [1:0] op;        // ALU operation code

    // DUT outputs (not randomized)
    logic [3:0] out;            // ALU result
    logic       c;              // Carry-out flag

    // Constraint to ensure all ALU operations
    // (0, 1, 2, 3) occur with equal probability
    constraint operation_t {
      op dist {[0:3] := 25};
    }

  endclass : ALU_transaction

endpackage : ALU_transaction_pkg
