`ifndef ALU_ITEM_SV      // <--- 1. Check if this flag is undefined
`define ALU_ITEM_SV      // <--- 2. Define the flag
import alu_pack::*;  // Import the package containing the enum
import uvm_pkg::*;
`include "uvm_macros.svh"
class seq_item extends uvm_sequence_item;
`uvm_object_utils(seq_item)
  randc bit [3:0] a;
  randc bit [3:0] b;
  operation_t op; 
    bit       c;
    bit [3:0] out;
    function new(string name= "seq_item");
    super.new(name); 
    endfunction
    function string convert2string();
return $sformatf ("A =%d ,B=%d , Op=%d ,C=%h ,out=%d ",
  a, b, op, c, out);
    endfunction

endclass
`endif