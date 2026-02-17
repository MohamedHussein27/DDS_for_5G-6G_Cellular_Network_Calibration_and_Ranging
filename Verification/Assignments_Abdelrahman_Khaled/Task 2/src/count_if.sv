// Interface definition for FSM count signals
// Used to connect DUT, driver, and monitor
interface count_if (input bit clk);

  // Control and configuration signals
  logic         rst_n;        // Active-low async reset
  logic         start;        // Start pulse
  logic         flag;         // Stop condition flag
  logic [15:0]  wait_timer;   // Programmable wait cycles

  // Status outputs
  logic        busy;         // FSM running indicator
  logic [4:0]  count_value;  // Count result

  // -----------------------------
  // Modport for the DUT
  // -----------------------------
  modport DUT (
    input  clk, rst_n, start, flag, wait_timer,
    output busy, count_value
  );

  // -----------------------------
  // Modport for the testbench driver
  // -----------------------------
  modport tb (
    input  clk, busy, count_value,
    output rst_n, start, flag, wait_timer
  );

  // -----------------------------
  // Modport for the monitor
  // -----------------------------
  modport monitor (
    input clk, rst_n, start, flag, wait_timer, busy, count_value
  );

endinterface : count_if
