interface alu_if;
    bit        clk;
    logic [3:0] a;
    logic [3:0] b;
    logic [1:0] op;
    logic       c;
    logic [3:0] out;
    
    initial begin
      forever begin
         #10;
         clk = ~clk;
      end
   end

endinterface