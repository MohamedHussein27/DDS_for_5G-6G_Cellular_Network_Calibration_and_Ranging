`include "counter_item.sv"
class counter_coverage;
    counter_item item;

    covergroup cg;
        cp_start:  coverpoint item.start;
        cp_flag:  coverpoint item.flag;
        cp_wait_timer: coverpoint item.wait_timer;
        
    endgroup

    function new();
        cg = new();
    endfunction

    task sample(counter_item t);
        this.item = t;
        cg.sample();
    endtask
endclass