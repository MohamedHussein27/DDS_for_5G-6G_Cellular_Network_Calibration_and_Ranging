`include "ALU_ENV.sv"
import shared_pkg::*;
module alu_top;

    bit clk;

    // clock generation
    initial begin
        clk = 0;
        forever #10 clk = ~clk;
    end

    alu_if aluif(clk);

    ALU dut (
        .rst_n(aluif.rst_n),
        .clk (aluif.clk  ),
        .a (aluif.a  ),
        .b (aluif.b  ),
        .op(aluif.op ),
        .out(aluif.out),
        .c (aluif.c  )
    );

    alu_env env;

    initial begin
        env = new();
        env.set_vif(aluif); // pass the real interface to env
        env.connect_mails(); // connect mailboxes
        env.run();
        $display("error count = %0d, correct count = %0d", error_count_out, correct_count_out);
        $stop;
    end

endmodule
