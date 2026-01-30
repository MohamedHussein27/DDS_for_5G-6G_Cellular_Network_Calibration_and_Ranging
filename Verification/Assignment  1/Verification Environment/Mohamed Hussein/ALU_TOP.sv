module alu_top();
    // interface
    alu_if aluif ();
    // dut
    ALU dut (aluif);
    // test
    alu_driver drv (aluif);
    // monitor
    alu_monitor MONITOR (aluif);
endmodule