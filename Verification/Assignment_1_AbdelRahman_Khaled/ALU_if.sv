// Interface definition for ALU signals
// Used to connect DUT, testbench, and monitor
interface ALU_if ();

  // ALU operation selector
  logic [1:0] op;

  // ALU input operands
  logic [3:0] a;
  logic [3:0] b;

  // ALU outputs
  logic       c;     // Carry-out flag
  logic [3:0] out;   // ALU result

  // Modport for the Design Under Test (DUT)
  // DUT reads inputs and drives outputs
  modport DUT (
    input  op, a, b,
    output c, out
  );

  // Modport for the testbench driver
  // Testbench drives inputs and reads outputs
  modport tb (
    output op, a, b,
    input  c, out
  );

  // Modport for the monitor
  // Monitor passively observes all signals
  modport monitor (
    input op, a, b, c, out
  );

endinterface : ALU_if
