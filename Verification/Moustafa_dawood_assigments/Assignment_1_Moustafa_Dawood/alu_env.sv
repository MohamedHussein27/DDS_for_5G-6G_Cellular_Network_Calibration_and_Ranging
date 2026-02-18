
import alu_pkg::*;

`include "alu_sequnce_item.sv"
`include "alu_sequnce.sv"
`include "alu_driver.sv"
`include "alu_monitor.sv"
`include "alu_scoreboard.sv"
`include "alu_coverage.sv"


class alu_env;

  alu_sequence seq;
  alu_driver     drv;
  alu_monitor    mon;
  alu_scoreboard sb;
  alu_coverage   cov;

  mailbox #(alu_seq_item) s2d;
  mailbox #(alu_seq_item) d2m;
  mailbox #(alu_seq_item) m2s;

  function new(virtual alu_bfm bfm);
    s2d = new();
    d2m = new();
    m2s = new();

    seq = new(s2d);
    drv = new(bfm, s2d, d2m);
    mon = new(d2m, m2s);
    sb  = new(m2s);
    cov = new(bfm);
  endfunction

  task run();
    fork
      seq.run();
      drv.run();
      mon.run();
      sb.run();
      cov.run();
    join_none
  endtask

endclass
