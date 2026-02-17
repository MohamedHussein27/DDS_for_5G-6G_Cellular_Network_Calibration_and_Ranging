`include "monitor.sv"
`include "driver.sv"
`include "sequencer.sv"
`include "alu_config.sv"
class agent extends uvm_agent;
`uvm_component_utils (agent)
sequencer sequencer_;
driver driver_;
monitor monitor_;
alu_config cfg;
uvm_analysis_port #(seq_item) agent_port;
function new(string name = "agent", uvm_component parent = null);
super.new(name, parent);
endfunction
function void build_phase(uvm_phase phase);
super.build_phase(phase);
if (!uvm_config_db#(alu_config)::get(this, "", "vif", cfg)) 
    `uvm_fatal("NOVIF", "Virtual interface not found in config DB")
sequencer_ = sequencer::type_id::create("sequencer_", this);
driver_ = driver::type_id::create("driver_", this);
monitor_ = monitor::type_id::create("monitor_", this);
agent_port = new("agent_port", this);
endfunction
function void connect_phase(uvm_phase phase);
super.connect_phase(phase);
driver_.vif = cfg.vif;
monitor_.vif = cfg.vif;
driver_.seq_item_port.connect(sequencer_.seq_item_export);
monitor_.mon_port.connect(agent_port);
endfunction

endclass