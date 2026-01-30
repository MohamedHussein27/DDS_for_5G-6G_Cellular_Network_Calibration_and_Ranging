interface alu_if();

    logic [3:0] a;
    logic [3:0] b;
    logic [1:0] op;
    logic       c;
    logic [3:0] out;

    modport DUT (
        input  a, b, op, c,
        output out
    );

    modport TEST (
        output a, b, op, c,
        input  out
    );

    modport MONITOR (
        input a, b, op, c, out
    );

endinterface