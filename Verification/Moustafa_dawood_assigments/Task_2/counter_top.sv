import counter_env_pkg::*;
`include "counter_if.sv"   // counter_bfm interface
`timescale 1ns/1ps
module top;

  // -----------------------------------
  // BFM
  // -----------------------------------
  counter_bfm bfm();

  logic clk;

  // -----------------------------------
  // Clock generation (TOP responsibility)
  // -----------------------------------
  initial begin
    clk = 0;
    forever #5 clk = ~clk;
  end

  assign bfm.clk = clk;

  // -----------------------------------
  // DUT
  // -----------------------------------
  count_fsm dut (
    .clk         (clk),
    .rst_n       (bfm.rst_n),
    .start       (bfm.start),
    .flag        (bfm.flag),
    .wait_timer  (bfm.wait_timer),
    .busy        (bfm.busy),
    .count_value (bfm.count_value)
  );
// assertions
  counter_assertions dut_asserts (.bfm(bfm));

  // -----------------------------------
  // Environment
  // -----------------------------------
  counter_env env;

  // -----------------------------------
  // Test flow
  // -----------------------------------
  initial begin
    env = new();          // build
    env.connect(bfm);    // connect
    env.run();           // run (build/connect/run like UVM)

    #20000;
    $stop;
  end

endmodule
