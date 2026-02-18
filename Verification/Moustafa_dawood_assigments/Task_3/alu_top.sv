import uvm_pkg::*;
import alu_test_pkg::*;
`include "uvm_macros.svh"

module alu_top ();

bit clk;

	initial begin
        clk = 0;             // Drive local 'clk'
        forever #1 clk = ~clk; 
    end

    alu_if alu_if_ins(clk); // making instance of the real interface
    ALU dut ( .clk   (alu_if_ins.clk),
    .rst_n (alu_if_ins.rst_n),
    .a     (alu_if_ins.operand_a),
    .b     (alu_if_ins.operand_b),
    .op    (alu_if_ins.op),
    .c     (alu_if_ins.c),
    .out   (alu_if_ins.out)); // we must go to dut and connect the real interface to the dut
    
    initial begin
		uvm_config_db#(virtual alu_if)::set(null, "uvm_test_top", "ALU_IF", alu_if_ins);
        //data type called virtual interface. This is a variable that can point to a static real interface instance

		//run_test("alu_test");
        run_test();
	end
endmodule