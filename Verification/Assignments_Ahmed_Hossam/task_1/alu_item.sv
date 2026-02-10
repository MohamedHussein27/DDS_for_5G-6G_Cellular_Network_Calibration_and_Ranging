`ifndef ALU_ITEM_SV      // <--- 1. Check if this flag is undefined
`define ALU_ITEM_SV      // <--- 2. Define the flag
import alu_pack::*;  // Import the package containing the enum
class alu_item;
  randc bit [3:0] a;
  randc bit [3:0] b;
  randc operation_t op; 
    bit       c;
    bit [3:0] out;
 constraint valid_op_c {
    op inside {
      alu_pack::ADD_op,
      alu_pack::XOR_op,
      alu_pack::AND_op,
      alu_pack::OR_op
    };
 }
endclass

`endif                  