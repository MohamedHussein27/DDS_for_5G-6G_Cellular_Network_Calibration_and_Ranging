interface seq_intf(input logic clk);
    // Signals
    logic rst_n;
    logic data_in;
    logic seq_detected; // Output from DUT, Input to TB

    modport DUT (
        input  clk, rst_n, data_in,
        output seq_detected
    );

    modport TEST (
        input  clk, seq_detected,
        output rst_n, data_in
    );
endinterface