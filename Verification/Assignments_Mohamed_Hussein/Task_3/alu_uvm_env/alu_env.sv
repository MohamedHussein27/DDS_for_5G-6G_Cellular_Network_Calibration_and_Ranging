package alu_env_pkg;
    import alu_config_obj_pkg::*;
    import alu_seq_item_pkg::*;
    import alu_sequencer_pkg::*;
    import alu_agent_pkg::*;
    import alu_subscriber_pkg::*;
    import alu_scoreboard_pkg::*;
    import uvm_pkg::*;
    `include "uvm_macros.svh"
    class alu_env extends uvm_env;
        `uvm_component_utils (alu_env)

        // agent, scoreboard and sub
        alu_agent agt;
        alu_scoreboard sb;
        alu_subscriber sub;

        // configuration object 
        alu_config_obj alu_cfg; // configuration object     
    
        // construction
        function new (string name = "alu_env", uvm_component parent = null);
            super.new(name, parent);
        endfunction


        function void build_phase (uvm_phase phase);
            super.build_phase(phase);
            agt = alu_agent::type_id::create("agt",this);            
            sb = alu_scoreboard::type_id::create("sb", this);
            sub = alu_subscriber::type_id::create("sub", this);
        endfunction

        // connection between agent and scoreboard and between agent and sub 
        function void connect_phase (uvm_phase phase);
            agt.agt_ap.connect(sb.sb_export);
            agt.agt_ap.connect(sub.sub_export);
        endfunction
    endclass
endpackage