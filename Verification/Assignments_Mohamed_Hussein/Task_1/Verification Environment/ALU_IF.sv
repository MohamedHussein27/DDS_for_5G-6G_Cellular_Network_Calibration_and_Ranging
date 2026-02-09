interface alu_if(clk);

    input clk;

    logic rst_n;
    logic [3:0] a;
    logic [3:0] b;
    logic [1:0] op;
    logic       c;
    logic [3:0] out;

    modport DUT (
        input clk, rst_n, a, b, op, c,
        output out
    );

    modport TEST (
        output rst_n, a, b, op, c,
        input clk, out
    );

    modport MONITOR (
        input clk, rst_n, a, b, op, c, out
    );

    // clocking block for monitor to sample data slightly after posedge of clock
    clocking mon_cb @(posedge clk);
        default input #1step output #1step; 
        output rst_n, a, b, op; 
        input out, c;
    endclocking

endinterface