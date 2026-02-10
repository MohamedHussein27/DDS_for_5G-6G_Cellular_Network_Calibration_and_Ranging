interface counter_if(clk);

    input clk;

    logic rst_n;
    logic start;
    logic [15:0] wait_timer;
    logic flag;
    logic busy;
    logic [4:0] count_value;

    // clocking block for monitor to sample data slightly after posedge of clock
    clocking mon_cb @(posedge clk);
        default input #1step output #1step; 
        output rst_n, start, wait_timer, flag; 
        input busy, count_value;
    endclocking

endinterface