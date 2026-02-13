interface alu_if(clk);

    input clk;

    logic rst_n;
    logic [3:0] a;
    logic [3:0] b;
    logic [1:0] op;
    logic       c;
    logic [3:0] out;

endinterface