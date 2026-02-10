`include "counter_env.sv" // This pulls in everything else!

module tb_top;
    
    // 1. Interface
    counter_if intf();
    
    // 2. DUT
    count_fsm dut (
        .clk(intf.clk),
        .rst_n(intf.rst_n),
        .start(intf.start),
        .wait_timer(intf.wait_timer),
        .flag(intf.flag),
        .busy(intf.busy),
        .count_value(intf.count_value)
    );

    // 3. Environment
    counter_env env;

    initial begin
        env = new();
            env.connecting(intf); // Connect env to interface
        
        $display("Starting Simulation...");
        env.run(); 
        #6000
    
        $finish;
    end
endmodule