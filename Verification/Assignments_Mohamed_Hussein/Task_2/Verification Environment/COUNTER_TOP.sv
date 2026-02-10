`include "COUNTER_ENV.sv"
import shared_pkg::*;
module counter_top;

    bit clk;

    // clock generation
    initial begin
        clk = 0;
        forever #10 clk = ~clk;
    end

    counter_if counterif(clk);

    count_fsm dut (
        .rst_n      (counterif.rst_n),
        .clk        (counterif.clk),
        .start      (counterif.start),
        .wait_timer (counterif.wait_timer),
        .flag       (counterif.flag ),
        .busy       (counterif.busy),
        .count_value(counterif.count_value)
    );
    counter_env env;

    initial begin
        env = new();
        env.counter_vif = counterif; // assigning real interface to env virtual interface
        env.connect(); // function used to connect mailboxes and virtual interfaces
        env.run();
        $stop;
    end

endmodule