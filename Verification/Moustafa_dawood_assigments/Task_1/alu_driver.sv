
/*
class alu_driver;
  

  virtual alu_bfm bfm;
  mailbox #(alu_seq_item) seq2drv_mb;
  mailbox #(alu_seq_item) drv2mon_mb;

  function new(
    virtual alu_bfm b,
    mailbox #(alu_seq_item) s2d,
    mailbox #(alu_seq_item) d2m
  );
    bfm = b;
    seq2drv_mb = s2d;
    drv2mon_mb = d2m;
  endfunction

  task run();
    alu_seq_item item;

    forever begin
      seq2drv_mb.get(item);

      bfm.operand_a = item.a;
      bfm.operand_b = item.b;
      bfm.op_set    = item.op;

      #1; // combinational settle
      item.result = bfm.result;

      drv2mon_mb.put(item);
    end
  endtask

endclass

*/
/*
class alu_driver;

  virtual alu_bfm bfm;
  mailbox #(alu_seq_item) seq2drv_mb;
  mailbox #(alu_seq_item) drv2mon_mb;

  function new(
    virtual alu_bfm b,
    mailbox #(alu_seq_item) s2d,
    mailbox #(alu_seq_item) d2m
  );
    bfm = b;
    seq2drv_mb = s2d;
    drv2mon_mb = d2m;
  endfunction

    // Task to drive reset 
  task drive_reset(bit rst_value = 0, int cycles = 2);
    // assert reset
    bfm.rst_n <= rst_value;
    repeat(cycles) @(negedge bfm.clk);
    // deassert reset
    bfm.rst_n <= 1'b1;
    @(negedge bfm.clk);
  endtask
// my main main run 
  task run();
    alu_seq_item item;

    forever begin
      seq2drv_mb.get(item);

      // Drive on negedge edge
      @(negedge bfm.clk);
      bfm.operand_a <= item.a;
      bfm.operand_b <= item.b;
      bfm.op_set    <= item.op;

      // Pass item forward 
      drv2mon_mb.put(item);
    end
  endtask

endclass
*/
/*
class alu_driver;

  virtual alu_bfm bfm;
  mailbox #(alu_seq_item) seq2drv_mb;
  mailbox #(alu_seq_item) drv2mon_mb;

  function new();
  endfunction

  // Manual connection of BFM and mailboxes
  function void connect(virtual alu_bfm b, mailbox #(alu_seq_item) s2d, mailbox #(alu_seq_item) d2m);
    bfm        = b;
    seq2drv_mb = s2d;
    drv2mon_mb = d2m;
  endfunction

  // Task to drive reset
  task drive_reset(logic rst_val, int cycles);
    bfm.rst_n <= rst_val;
    repeat(cycles) @(posedge bfm.clk);
    bfm.rst_n <= 1'b1;  // release reset
    @(posedge bfm.clk);
  endtask

  task run();
    alu_seq_item item;

    forever begin
      seq2drv_mb.get(item);

      // Drive DUT inputs on negedge (typical for combinational DUT)
      @(negedge bfm.clk);
      bfm.operand_a <= item.a;
      bfm.operand_b <= item.b;
      bfm.op_set    <= item.op;

      // Put item for monitor
      drv2mon_mb.put(item);
    end
  endtask

endclass
*/


//import alu_pkg::*;
//`include "alu_seq_item.sv"
package alu_driver_pkg;

import alu_sequnce_item_pkg::*;
class alu_driver;

  virtual alu_bfm bfm;
  mailbox #(alu_seq_item) seq2drv_mb;
  mailbox #(alu_seq_item) drv2mon_mb;

  function new();
  endfunction

  //  Reset task with NO argument
  task drive_reset(); 
    $display("DRIVER Asserting Reset");
    bfm.rst_n <= 0;                // Drive Reset Low
    repeat(2) @(posedge bfm.clk);  // Wait 2 Clock Cycles
    bfm.rst_n <= 1;                // Release Reset
    @(posedge bfm.clk);            // Wait one more cycle
    $display("DRIVER Reset Done");
  endtask
  

  task run();
    alu_seq_item item;

    forever begin
      seq2drv_mb.get(item);

      // Drive at negedge not puting anything 
      @(negedge bfm.clk);
      bfm.operand_a <= item.a;
      bfm.operand_b <= item.b;
      bfm.rst_n     <= item.rst_n;
      bfm.op_set    <= item.op;

      // Send item to monitor
     // drv2mon_mb.put(item);
    end
  endtask

endclass

endpackage