
interface counter_bfm;

  // Clock & Reset
  logic        clk;
  logic        rst_n;

  // DUT inputs
  logic        start;
  logic [15:0] wait_timer;
  logic        flag;

  // DUT outputs
  logic        busy;
  logic [4:0]  count_value;

endinterface
