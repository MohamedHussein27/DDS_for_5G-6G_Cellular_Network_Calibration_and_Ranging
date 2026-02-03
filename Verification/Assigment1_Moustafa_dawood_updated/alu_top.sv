//`include "alu_env.sv"
/*
module top;
  import alu_pkg::*;

  alu_bfm bfm();
  logic clk;
  logic rst_n;

  // clk generation must done in top 
  initial begin
    clk=0;
    forever #5 clk = ~clk;
  end


  ALU dut (
    .clk (clk),
    .rst_n (bfm.rst_n),
    .a   (bfm.operand_a),
    .b   (bfm.operand_b),
    .op  (bfm.op_set),
    .c   (bfm.c),
    .out (bfm.out)
  );

  alu_env env;

  initial begin
 env = new();    // just build objects
env.connect(bfm);       // manually connect all
env.run();              // start all 
//(we here simulating build and connect and run phase of uvm )

    #20000
     $stop;
  end

endmodule
*/
// ALU project all includes

import alu_env_pkg::*;
`include "alu_if.sv" // Interface file



module top;


  // BFM
  alu_bfm bfm();
  logic clk;
  // Clock generation
  initial begin
    clk = 0;
    forever #5 clk = ~clk;
  end
assign bfm.clk = clk;
  // DUT
  ALU dut (
    .clk (clk),
    .rst_n (bfm.rst_n),
    .a   (bfm.operand_a),
    .b   (bfm.operand_b),
    .op  (bfm.op_set),
    .c   (bfm.c),
    .out (bfm.out)
  );

  // Environment
  alu_env env;

  

  initial begin
    env = new();
    env.connect(bfm);   // connect manually
    env.run();          // run all components
    #20000;
   $stop;
  end

endmodule



/*
alu_pkg.sv
alu_if.sv

alu_seq_item.sv
alu_sequence.sv
alu_driver.sv
alu_monitor.sv
alu_scoreboard.sv
alu_coverage.sv
alu_env.sv

top.sv
*/
