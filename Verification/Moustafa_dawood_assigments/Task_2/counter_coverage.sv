
package counter_coverage_pkg;
import counter_seq_item_pkg::*;

class counter_coverage;

  mailbox #(counter_seq_item) m2cov_mb;
  counter_seq_item item;

  // Local sampled variables
  bit start;
  bit flag;
  bit busy;
  bit [15:0] wait_timer;
  bit [4:0] count_value;

  covergroup counter_cg;
    cp_start : coverpoint start;
    cp_flag  : coverpoint flag;
    cp_busy  : coverpoint busy;

coverpoint wait_timer {
    bins small_bin = {[0:10]};
    bins mid_bin   = {[11:100]};
    bins large_bin = {[101:$]};
}


    cp_count : coverpoint count_value {
      bins zero = {0};
      bins mid  = {[1:15]};
      bins max  = {[16:31]};
    }

    //  Cross coverages
    cross cp_start, cp_busy;
    cross cp_flag, cp_busy;
  endgroup

  function new();
    counter_cg = new();
  endfunction

  task run();
    forever begin
      m2cov_mb.get(item);

      start       = item.start;
      flag        = item.flag;
      busy        = item.busy;
      wait_timer  = item.wait_timer;
      count_value = item.count_value;

      counter_cg.sample();
    end
  endtask

endclass
endpackage
