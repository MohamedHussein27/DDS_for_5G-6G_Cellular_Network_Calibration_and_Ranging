`include "alu_env.sv" // This pulls in everything else!

module tb_top;
    
    // 1. Interface
    alu_if intf();
    
    // 2. DUT
    ALU dut (
        .clk(intf.clk),
        .a(intf.a),
        .b(intf.b),
        .op(intf.op),
        .c(intf.c),
        .out(intf.out)
    );

    // 3. Environment
    alu_env env;

    initial begin
        env = new(intf);
        
        $display("Starting Simulation...");
        env.run(); // Run 3000 tests
        #20500
    
        $finish;
    end
endmodule