// Include the Class Definitions we wrote previously
`include "Seq_Detector_Scoreboard.sv" 

module tb_top;

    // Clock Generation
    bit clk;
    // clock generation
    initial begin
        clk = 0;
        forever #5 clk = ~clk;
    end

    // Instantiate Interface
    seq_intf intf(clk);

    // Instantiate DUT (Connect to Interface)
    seq_detector DUT (
        .clk(intf.clk),
        .rst_n(intf.rst_n),
        .data_in(intf.data_in),
        .seq_detected(intf.seq_detected)
    );

    // Declare Scoreboard Handle
    my_scoreboard sb;

    // Main Test Process
    initial begin
        intf.rst_n = 0;
        intf.data_in = 0;

        sb = new();
        sb.set_vif(intf);

        fork
            sb.run();
        join_none

        #20;
        intf.rst_n = 1;
        $display("--- Reset Released ---");

        
        // Sequence: 1 -> 0 -> 1 (Detect)
        @(negedge clk) intf.data_in = 1; 
        @(negedge clk) intf.data_in = 0; 
        @(negedge clk) intf.data_in = 1; // Expect "Match: OK" for S101
        
        // Sequence: 1 (Overlap Detect) 
        @(negedge clk) intf.data_in = 1; // Expect "Match: OK" for S1 (from overlap)
        @(negedge clk) intf.data_in = 0; 
        //@(negedge clk) intf.data_in = 0; 

        intf.rst_n = 0;
        intf.data_in = 0;
        #20;
        intf.rst_n = 1;
        $display("--- Reset Released again ---");

        // Sequence: 1 -> 0 -> 1 (Detect)
        @(negedge clk) intf.data_in = 1; 
        @(negedge clk) intf.data_in = 0; 
        @(negedge clk) intf.data_in = 1; // Expect "Match: OK" for S101

        // Wait a bit and finish
        #50;
        $display("--- Test Completed ---");
        $finish;
    end


endmodule