`timescale 1ns/1ps
package counter_monitor_pkg;
import counter_seq_item_pkg::*;
class counter_monitor;
  
  virtual counter_bfm bfm;
  mailbox #(counter_seq_item) mon2sb_mb;
  mailbox #(counter_seq_item) m2cov_mb;
   
  
    function new();
        
    endfunction //new()

    task run();
        counter_seq_item item;
        forever begin
            //drv2mon_mb.get(item);
            item = new();
            @(posedge bfm.clk);
            #0.00000001; // small delay to ensure stable values
            item.rst_n = bfm.rst_n;
            item.start = bfm.start;
            item.flag = bfm.flag;
            item.wait_timer = bfm.wait_timer;
            item.busy = bfm.busy;
            item.count_value = bfm.count_value;
            mon2sb_mb.put(item);
            m2cov_mb.put(item);
        end

        endtask
endclass //className
    
endpackage