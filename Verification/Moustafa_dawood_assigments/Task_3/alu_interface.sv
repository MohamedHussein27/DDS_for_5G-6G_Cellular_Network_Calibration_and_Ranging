import alu_pkg::*;
interface alu_if(input logic clk);


    // ALU BFM signals
   
    logic [3:0]   operand_a;
    logic [3:0]   operand_b;
    logic [4:0]   result;
    opcode_t      op_set;
    // clk and reset signals added 
    //logic clk;
    logic rst_n;



       //  OUTPUTS TO CONNECT TO DUT:
    logic [1:0]   op;        // Connect to DUT op input
    logic         c;         // Connect from DUT c output
    logic [3:0]   out;       // Connect from DUT out output
    
    // Connect DUT outputs to result
    assign result = {c, out};
    assign op = op_set;


   
    
endinterface //interfacename