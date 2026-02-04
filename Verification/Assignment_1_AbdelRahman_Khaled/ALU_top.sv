import ALU_env_pkg::*; 

module ALU_top;

  // Clock signal
  bit clk;

  // 2. Instantiate Interface (This holds the actual signals a, b, op, etc.)
  ALU_if alu_if(clk);

  // 3. Instantiate DUT
  // CONNECT DIRECTLY TO THE INTERFACE SIGNALS
  ALU dut (
    .clk  (alu_if.clk),
    .a   (alu_if.a),    
    .b   (alu_if.b),
    .op  (alu_if.op),
    .out (alu_if.out),
    .c   (alu_if.c)
  );

  // Environment handle
  ALU_env env;

  // Clock generation
  initial begin
    clk = 0;
    forever #1 clk = ~clk;  // 100MHz clock
  end

  // Test Flow
  initial begin
    // Create environment and pass virtual interfaces
    env = new(alu_if.tb, alu_if.monitor);
    
    env.build();   // Build hierarchy
    env.run();     // Start simulation
    $stop;
  end

endmodule