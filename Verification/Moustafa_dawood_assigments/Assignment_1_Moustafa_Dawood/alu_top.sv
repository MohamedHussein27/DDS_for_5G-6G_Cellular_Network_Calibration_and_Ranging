`include "alu_env.sv"
module top;
  import alu_pkg::*;

  alu_bfm bfm();

  ALU dut (
    .a   (bfm.operand_a),
    .b   (bfm.operand_b),
    .op  (bfm.op_set),
    .c   (bfm.c),
    .out (bfm.out)
  );

  alu_env env;

  initial begin
    env = new(bfm);
    env.run();
    #20000
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
