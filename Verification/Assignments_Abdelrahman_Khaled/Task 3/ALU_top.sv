`timescale 1ns/10ps

import uvm_pkg::*;
`include "uvm_macros.svh"

import ALU_test_pkg::*;

module ALU_top;

  //--------------------------------------------------------------------------
  // Clock Generation
  //--------------------------------------------------------------------------
  bit clk;

  initial begin
    clk = 0;
    forever #5 clk = ~clk;   // 10ns clock period
  end

  //--------------------------------------------------------------------------
  // Interface Instance
  //--------------------------------------------------------------------------
  ALU_if vif(clk);

  //--------------------------------------------------------------------------
  // DUT Instance
  //--------------------------------------------------------------------------
  ALU dut (
    .clk (vif.clk),
    .a   (vif.a),
    .b   (vif.b),
    .op  (vif.op),
    .out (vif.out),
    .c   (vif.c)
  );

  //--------------------------------------------------------------------------
  // UVM Configuration + Test Start
  //--------------------------------------------------------------------------
  initial begin

    // Provide virtual interface to UVM test
    uvm_config_db#(virtual ALU_if)::set(
      null,
      "uvm_test_top",
      "ALU_IF",
      vif
    );

    // Run selected test
    run_test();   // Recommended (choose test from command line)

  end

endmodule
