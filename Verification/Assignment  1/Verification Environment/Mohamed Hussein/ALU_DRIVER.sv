import alu_seq_item_pkg::*;
import shared_pkg::*;
module alu_driver (alu_if.TEST aluif);
    logic [3:0] a;
    logic [3:0] b;
    logic [1:0] op;
    logic [3:0] out;
    logic c;
    // assigning signals to be interfaced
    // inputs
    assign out = aluif.out;
    assign c = aluif.c;
    // outputs
    assign aluif.a = a;
    assign aluif.b = b;
    assign aluif.op = op;

    // class object
    alu_seq_item seq_dr = new;

    // stimulus
    initial begin
        repeat(15000) begin
            assert(seq_dr.randomize());
            /*aluif.a = seq_dr.a;
            aluif.b = seq_dr.b;
            aluif.op = seq_dr.op;*/
            a = seq_dr.a;
            b = seq_dr.b;
            op = seq_dr.op;
            #1;
        end
        test_finished = 1; // end of test
    end
endmodule