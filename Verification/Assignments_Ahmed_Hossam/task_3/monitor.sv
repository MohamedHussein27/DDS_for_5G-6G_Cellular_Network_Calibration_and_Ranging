`include "seq_item.sv"
class monitor extends uvm_monitor;
`uvm_component_utils (monitor)
virtual alu_if vif;
seq_item item;
uvm_analysis_port #(seq_item) mon_port;
function new(string name = "monitor", uvm_component parent = null);
super.new(name, parent);
endfunction
function void build_phase(uvm_phase phase);
super.build_phase(phase);
mon_port = new("mon_port", this);
endfunction
task run_phase(uvm_phase phase);
super.run_phase(phase);
forever begin
    item =seq_item::type_id::create("item");
    @(negedge vif.clk);
    item.a  = vif.a;
    item.b  = vif.b;
    void'($cast(item.op, vif.op)); // Cast 2-bit wire to enum
    item.out = vif.out;
    item.c   = vif.c;
    mon_port.write(item);
    `uvm_info("run_phase",item.convert2string(),UVM_HIGH)
end
endtask
endclass