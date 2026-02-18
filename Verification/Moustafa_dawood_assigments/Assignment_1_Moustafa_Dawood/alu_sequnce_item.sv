


class alu_seq_item;
  
  randc bit [3:0] a;
  randc bit [3:0] b;
  randc opcode_t  op;  // 

  bit [4:0] result;   // Actual DUT output

    //  constraint
 constraint valid_op_c {
    op inside {
      alu_pkg::ADD_op,
      alu_pkg::XOR_op,
      alu_pkg::AND_op,
      alu_pkg::OR_op
    };
 }
endclass
