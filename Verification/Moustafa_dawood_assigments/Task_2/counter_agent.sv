//import alu_pkg::*;
package counter_agent_pkg;
import counter_sequence_pkg::*;
import counter_driver_pkg::*;
import counter_monitor_pkg::*;
import counter_seq_item_pkg::*;
class counter_agent;

  // Components
  counter_sequence seq;
  counter_driver   drv;
  counter_monitor  mon;

  // Mailboxes
  mailbox #(counter_seq_item) s2d;
 // mailbox #(counter_seq_item) d2m;
  mailbox #(counter_seq_item) m2cov_mb;
  mailbox #(counter_seq_item) m2sb;


  // Constructor: just create objects and mailboxes
  function new();
    s2d = new();
   
    m2sb= new(); 
    m2cov_mb = new();
    seq = new();
    drv = new();
    mon = new();
  endfunction

  // Manual connections
  function void connect(virtual counter_bfm bfm);
    // Connect sequence to driver
    seq.connect(s2d);

    // Connect driver
    drv.bfm        = bfm;
    drv.seq2drv_mb = s2d;
   

    // Connect monitor
    mon.bfm        = bfm;
  
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