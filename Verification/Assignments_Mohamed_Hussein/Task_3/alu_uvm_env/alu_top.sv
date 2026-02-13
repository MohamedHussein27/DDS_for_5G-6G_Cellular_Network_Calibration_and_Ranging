import alu_test_pkg::*;
import alu_env_pkg::*;
import uvm_pkg::*;
`include "uvm_macros.svh"

module alu_top();
    bit clk;

    initial begin
        forever #1 clk = ~clk;
    end

    alu_if aluif (clk);
    ALU dut (
        .rst_n(aluif.rst_n),
        .clk (aluif.clk  ),
        .a (aluif.a  ),
        .b (aluif.b  ),
        .op(aluif.op ),
        .out(aluif.out),
        .c (aluif.c  )
    );

    // setting the virtual interface to be accessible by the test
    initial begin
        uvm_config_db #(virtual alu_if)::set(null, "uvm_test_top", "alu_V", aluif); // main interface
        run_test ("alu_test");
    end
endmodule