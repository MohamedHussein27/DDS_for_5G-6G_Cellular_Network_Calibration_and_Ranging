// Interface definition for ALU signals
// Used to connect DUT, driver, and monitor
interface ALU_if (input bit clk);

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
    input  clk,
    input  op, a, b,
    output c, out
  );
endinterface : ALU_if