import ALU_tb_pkg::*;
import ALU_monitor_pkg::*;

module ALU_top();

  // Instantiate the ALU interface
  ALU_if ALUif();

  // Instantiate the ALU DUT, connecting it to the interface via the DUT modport
    ALU DUT (
        .a(ALUif.a),
        .b(ALUif.b),
        .op(ALUif.op),
        .out(ALUif.out),
        .c(ALUif.c)
    );

  // Instantiate the ALU testbench, connecting it to the interface via the TB modport
  ALU_tb tb;

  initial begin
    tb = new(ALUif.tb);  // pass virtual interface
    tb.run();
  end


  // Instantiate the ALU monitor, connecting it to the interface via the monitor modport
  ALU_monitor mon;

  initial begin
    mon = new(ALUif.monitor); // pass virtual interface
    mon.run();
  end


endmodule : ALU_top