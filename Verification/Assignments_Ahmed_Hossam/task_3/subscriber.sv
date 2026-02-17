`include "seq_item.sv"
class subscriber extends uvm_component;
`uvm_component_utils(subscriber)
uvm_analysis_export #(seq_item) cov_export;
uvm_tlm_analysis_fifo #(seq_item) cov_fifo;
seq_item item;
covergroup cg;
        cp_a:  coverpoint item.a;
        cp_b:  coverpoint item.b;
        cp_op: coverpoint item.op;
        cross_all: cross cp_a, cp_b, cp_op;
    endgroup
function new(string name = "subscriber", uvm_component parent = null);
super.new(name, parent);
cg=new();
endfunction
function void build_phase(uvm_phase phase);
super.build_phase(phase);
cov_export = new("cov_export", this);
cov_fifo = new("cov_fifo", this);
endfunction
function void connect_phase(uvm_phase phase);
super.connect_phase(phase);
cov_export.connect(cov_fifo.analysis_export);
endfunction
task run_phase(uvm_phase phase);
super.run_phase(phase);
forever begin
    cov_fifo.get(item);
    cg.sample();
    `uvm_info("run_phase",item.convert2string(),UVM_HIGH)
end
endtask
endclass