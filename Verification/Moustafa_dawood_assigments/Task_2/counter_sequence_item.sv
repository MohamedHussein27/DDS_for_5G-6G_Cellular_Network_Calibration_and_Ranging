package counter_seq_item_pkg;


class counter_seq_item;

  // Reset
  randc bit        rst_n;

  // DUT inputs
  randc bit        start;
  randc bit        flag;
  randc bit [15:0] wait_timer;

  // DUT outputs (sampled by monitor)
  logic           busy;
  logic [4:0]     count_value;

  
  // Constraints
  

  // Reset rarely asserted
  constraint rst_c {
    rst_n dist {0 := 1, 1 := 100};
  }

  //  avoid zero wait_timer
  constraint wait_timer_c {
    wait_timer inside {[1:16'hFFFF]};
  }

endclass

endpackage
