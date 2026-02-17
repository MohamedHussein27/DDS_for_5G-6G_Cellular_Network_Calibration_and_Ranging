import count_env_pkg::*;

`timescale 1ns / 1ps

module count_top;

  // -----------------------------
  // Clock
  // -----------------------------
  bit clk;

  // -----------------------------
  // Interface instantiation
  // -----------------------------
  count_if count_if (clk);

  // -----------------------------
  // DUT instantiation
  // -----------------------------
  count_fsm dut (
    .clk         (count_if.clk),
    .rst_n       (count_if.rst_n),
    .start       (count_if.start),
    .flag        (count_if.flag),
    .wait_timer  (count_if.wait_timer),
    .busy        (count_if.busy),
    .count_value (count_if.count_value)
  );

  // -----------------------------
  // Environment handle
  // -----------------------------
  count_env env;

  // -----------------------------
  // Clock generation
  // -----------------------------
  initial begin
    clk = 0;
    forever #2 clk = ~clk;   // 4 time-unit clock period
  end

  // -----------------------------
  // Test flow
  // -----------------------------
  initial begin
    // Create environment (empty constructor)
    env = new();

    // Connect environment to interfaces
    env.connect(count_if.tb, count_if.monitor);

    // Run test
    env.run();

    $stop;
  end

endmodule : count_top
