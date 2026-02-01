interface alu_bfm;

import alu_pkg::*;

    // ALU BFM signals
   
    logic [3:0]   operand_a;
    logic [3:0]   operand_b;
    logic [4:0]   result;
    opcode_t      op_set;


       //  OUTPUTS TO CONNECT TO DUT:
    logic [1:0]   op;        // Connect to DUT op input
    logic         c;         // Connect from DUT c output
    logic [3:0]   out;       // Connect from DUT out output
    
    // Connect DUT outputs to result
    assign result = {c, out};
    assign op = op_set;

   

   
/*
    task  send_op(input logic [3:0 ] ia,input logic  [3:0 ] ib, input opcode_t iop,output  logic [4:0] alu_result);
       @(negedge clk);
       operand_a = ia;
       operand_b = ib;
       op_set    = iop;
       alu_result = result;
       //@(posedge clk);
        
    endtask 
  */ 

    
endinterface //alu_bfm



