
`include "test.sv"
module top ();
alu_if alu_if_();
ALU dut (
    .clk(alu_if_.clk),
    .a(alu_if_.a),
    .b(alu_if_.b),
    .op(alu_if_.op),
    .c(alu_if_.c),
    .out(alu_if_.out)
);

initial begin
   uvm_config_db#(virtual alu_if)::set(null, "uvm_test_top", "ALU_IF", alu_if_); 
   run_test("");
end
endmodule