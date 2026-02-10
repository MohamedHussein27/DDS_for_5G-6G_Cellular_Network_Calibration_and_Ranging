interface counter_if;
    bit                    rst_n;
    bit                    clk;
    logic                  start;
    logic        [15:0]    wait_timer;
    logic                  flag;
    logic                  busy;
    logic        [4:0]     count_value; 

    initial begin
      forever begin
         #10;
         clk = ~clk;
      end
   end

endinterface