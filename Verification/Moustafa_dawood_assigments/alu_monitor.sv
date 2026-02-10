/*
class alu_monitor;

  mailbox #(alu_seq_item) drv2mon_mb;
  mailbox #(alu_seq_item) mon2sb_mb;

  function new(
    mailbox #(alu_seq_item) d2m,
    mailbox #(alu_seq_item) m2s
  );
    drv2mon_mb = d2m;
    mon2sb_mb  = m2s;
  endfunction

  task run();
    alu_seq_item item;
    forever begin
      drv2mon_mb.get(item);
      mon2sb_mb.put(item);
    end
  endtask

endclass
*/
//import alu_pkg::*;
//`include "alu_seq_item.sv"
package alu_monitor_pkg;
  
import alu_sequnce_item_pkg::*;
class alu_monitor;

  mailbox #(alu_seq_item) drv2mon_mb;
  mailbox #(alu_seq_item) mon2sb_mb;
  mailbox #(alu_seq_item) m2cov_mb;
  
  virtual alu_bfm bfm;

  function new();
  
  endfunction

  task run();
    alu_seq_item item;

    forever begin
      //drv2mon_mb.get(item);
    item = new();

      // Sample at posedge
      @(posedge bfm.clk);
      item.a     = bfm.operand_a;
      item.b     = bfm.operand_b;
      item.op    = bfm.op_set;
      item.rst_n = bfm.rst_n;
      item.result = {bfm.c, bfm.out};
      mon2sb_mb.put(item);
      m2cov_mb.put(item);

    end
  endtask

endclass
endpackage
