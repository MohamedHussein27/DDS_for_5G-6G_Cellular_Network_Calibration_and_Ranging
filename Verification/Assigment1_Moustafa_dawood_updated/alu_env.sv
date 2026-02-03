


/*
class alu_env;

  alu_sequence seq;
  alu_driver     drv;
  alu_monitor    mon;
  alu_scoreboard sb;
  alu_coverage   cov;

  mailbox #(alu_seq_item) s2d;
  mailbox #(alu_seq_item) d2m;
  mailbox #(alu_seq_item) m2s;

  // making the Constructor just create objects and mailboxes
  function new();
    s2d = new();
    d2m = new();
    m2s = new();

    seq = new();    // no mailbox 
    drv = new();    // no virtual interface
    mon = new();
    sb  = new();
    cov = new();
  endfunction

  // Manual connect procedure
  function void connect(virtual alu_bfm bfm);
    // Connect sequence to driver mailbox
    seq.seq2drv_mb = s2d;

    // Connect driver
    drv.bfm = bfm;
    drv.seq2drv_mb = s2d;
    drv.drv2mon_mb = d2m;

    // Connect monitor
    mon.drv2mon_mb = d2m;
    mon.mon2sb_mb  = m2s;

    // Connect scoreboard
    sb.mon2sb_mb = m2s;

    // Connect coverage
    cov.bfm = bfm;
  endfunction
// Runing sequnces here

task reset_sequence();
  drv.drive_reset(0, 5);  // drive reset low for 5 cycles
endtask

  // Run all components
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
*/
//import alu_pkg::*;
//`include "alu_agent.sv"
//`include "alu_scoreboard.sv"
//`include "alu_coverage.sv"

package alu_env_pkg;

import alu_agent_pkg::*;
import alu_scoreboard_pkg::*;
import alu_coverage_pkg::*;
import alu_sequnce_item_pkg::*;


class alu_env;

  // Components
  alu_agent      agent;
  alu_scoreboard sb;
  alu_coverage   cov;


  function new();
    agent = new();
    sb    = new();
    cov   = new();
  endfunction

  // Manual connections
  function void connect(virtual alu_bfm bfm);
    agent.connect(bfm);

    // Connect scoreboard
    sb.mon2sb_mb = agent.m2sb;
    


    // Connect coverage
    cov.bfm = bfm;
  endfunction

  // Run environment
  task run();
    fork
      agent.run();
      sb.run();
      cov.run();
    join_none
  endtask

endclass


endpackage