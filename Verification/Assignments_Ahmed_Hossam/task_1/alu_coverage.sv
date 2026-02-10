`include "alu_item.sv"
class alu_coverage;
    alu_item item;

    covergroup cg;
        cp_a:  coverpoint item.a;
        cp_b:  coverpoint item.b;
        cp_op: coverpoint item.op;
        cross_all: cross cp_a, cp_b, cp_op;
    endgroup

    function new();
        cg = new();
    endfunction

    task sample(alu_item t);
        this.item = t;
        cg.sample();
    endtask
endclass