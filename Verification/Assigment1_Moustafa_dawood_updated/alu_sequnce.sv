/*
class alu_sequence;

  mailbox #(alu_seq_item) seq2drv_mb;

  function new(mailbox #(alu_seq_item) mb);
    seq2drv_mb = mb;
  endfunction

  task run();
    alu_seq_item item;

    repeat (10000) begin
      item = new();
      assert(item.randomize());
      seq2drv_mb.put(item);
    end
  endtask

endclass
*/



/*
class alu_sequence;

  mailbox #(alu_seq_item) seq2drv_mb;

  function new(mailbox #(alu_seq_item) mb);
    seq2drv_mb = mb;
  endfunction

  task run();
    alu_seq_item item;

    for (int op = 0; op < 4; op++) begin
      for (int a = 0; a < 16; a++) begin
        for (int b = 0; b < 16; b++) begin
          item = new();
          item.op = opcode_t'(op);  
          item.a  = a;
          item.b  = b;
          seq2drv_mb.put(item);
        end
      end
    end
  endtask

endclass
*/
//import alu_pkg::*;
//`include "alu_seq_item.sv"


package alu_sequnce_pkg;
import alu_pkg::*;
import alu_sequnce_item_pkg::*;
class alu_sequence;

  mailbox #(alu_seq_item) seq2drv_mb;

  function new();
  endfunction

  // Connect mailbox manually
  function void connect(mailbox #(alu_seq_item) mb);
    seq2drv_mb = mb;
  endfunction

  // Task to generate all combinations
  task run();
    alu_seq_item item;
    for (int op = 0; op < 4; op++)
      for (int a = 0; a < 16; a++)
        for (int b = 0; b < 16; b++) begin
          item = new();
          item.op = opcode_t'(op);
          item.a  = a;
          item.b  = b;
          item.rst_n = 1; 
          seq2drv_mb.put(item);
        end
  endtask

  //  reset sequence task
  /*
  task rst_task();
    rst_n = 0;
    repeat (2) @(posedge clk);
    rst_n = 1;
    @(posedge clk);
  endtask
*/
endclass

endpackage
