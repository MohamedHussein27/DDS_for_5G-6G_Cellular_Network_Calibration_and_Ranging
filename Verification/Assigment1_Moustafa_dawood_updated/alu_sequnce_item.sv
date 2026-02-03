// i will not make changes here for adding clk and rst_n
package alu_sequnce_item_pkg;
import alu_pkg::*;
class alu_seq_item;
  

  rand bit [3:0] a;
  rand bit [3:0] b;
  rand opcode_t  op;
  rand bit        rst_n;
  bit [4:0] result;   // Actual DUT output

  // Constraint
  constraint valid_op_c {
    op inside { alu_pkg::ADD_op, alu_pkg::XOR_op, alu_pkg::AND_op, alu_pkg::OR_op };
  }
  constraint rst_c {
    rst_n dist {0 := 1, 1 := 100}; // 1/100 chance of reset being low
  }
endclass
endpackage