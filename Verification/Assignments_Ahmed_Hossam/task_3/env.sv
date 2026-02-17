`include "agent.sv"
`include "subscriber.sv"
`include "scoreboard.sv"
`include "seq_item.sv"
class env extends uvm_env; 
  `uvm_component_utils(env)
  agent agent_;
    subscriber subscriber_;
    scoreboard scb;
    function new(string name = "env", uvm_component parent = null);
        super.new(name, parent);
    endfunction
    function void build_phase(uvm_phase phase);
        super.build_phase(phase);
        agent_ = agent::type_id::create("agent_", this);
        subscriber_ = subscriber::type_id::create("subscriber_", this);
        scb = scoreboard::type_id::create("scb", this);
    endfunction
    function void connect_phase(uvm_phase phase);
    agent_.agent_port.connect(subscriber_.cov_export);
    agent_.agent_port.connect(scb.score_export);
    endfunction 
endclass