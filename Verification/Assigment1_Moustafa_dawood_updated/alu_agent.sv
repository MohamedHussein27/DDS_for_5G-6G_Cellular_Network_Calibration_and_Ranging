//import alu_pkg::*;
package alu_agent_pkg;
import alu_sequnce_pkg::*;
import alu_driver_pkg::*;
import alu_monitor_pkg::*;
import alu_sequnce_item_pkg::*;
class alu_agent;

  // Components
  alu_sequence seq;
  alu_driver   drv;
  alu_monitor  mon;

  // Mailboxes
  mailbox #(alu_seq_item) s2d;
 // mailbox #(alu_seq_item) d2m;
  mailbox #(alu_seq_item) m2cov_mb;
  mailbox #(alu_seq_item) m2sb;


  // Constructor: just create objects and mailboxes
  function new();
    s2d = new();
   // d2m = new();
       m2sb= new(); 
    seq = new();
    drv = new();
    mon = new();
  endfunction

  // Manual connections
  function void connect(virtual alu_bfm bfm);
    // Connect sequence to driver
    seq.connect(s2d);

    // Connect driver
    drv.bfm        = bfm;
    drv.seq2drv_mb = s2d;
   // drv.drv2mon_mb = d2m;

    // Connect monitor
    mon.bfm        = bfm;
  //  mon.drv2mon_mb = d2m;
    mon.m2cov_mb   = m2cov_mb;
    mon.mon2sb_mb  = m2sb;
  endfunction

  // Run agent: reset first, then normal sequence
  task run();
    // Reset DUT 
   drv.drive_reset();

    // Run sequence + driver + monitor
    fork
      seq.run();
      drv.run();
      mon.run();
    join_none
  endtask

endclass
endpackage