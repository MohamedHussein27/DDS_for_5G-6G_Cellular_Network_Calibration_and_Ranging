package alu_sequence_item_pkg;
import uvm_pkg::*;
import alu_pkg::*;

`include "uvm_macros.svh"

class alu_sequence_item extends uvm_sequence_item;
  `uvm_object_utils(alu_sequence_item)

  
    function new(string name = "alu_sequence_item");
      super.new(name);
        
    endfunction //new()

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


endclass //className extends superClass


    
endpackage