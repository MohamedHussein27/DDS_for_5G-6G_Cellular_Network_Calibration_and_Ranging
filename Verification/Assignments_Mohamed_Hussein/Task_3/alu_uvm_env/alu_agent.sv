package alu_agent_pkg;
    import alu_config_obj_pkg::*;
    import alu_driver_pkg::*;
    import alu_monitor_pkg::*;
    import alu_seq_item_pkg::*;
    import alu_sequencer_pkg::*;
    import uvm_pkg::*;
    `include "uvm_macros.svh"
    class alu_agent extends uvm_agent;
        `uvm_component_utils(alu_agent)

        alu_config_obj alu_cfg; // configuration object
        alu_driver drv; // driver
        alu_monitor mon; // monitor
        alu_sequencer sqr; // sequencer
        
        uvm_analysis_port #(alu_seq_item) agt_ap; // will be used to connect scoreboard and coverage collector

        function new (string name = "alu_agent", uvm_component parent = null);
            super.new(name, parent);
        endfunction

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);
            if(!uvm_config_db #(alu_config_obj)::get(this,"","CFG", alu_cfg)) begin
                    `uvm_fatal("build_phase","agent error");
            end
            // creation
            drv = alu_driver::type_id::create("driver", this);
            mon = alu_monitor::type_id::create("mon", this);
            sqr = alu_sequencer::type_id::create("sqr", this);
            agt_ap = new("agt_ap", this); // connection point
        endfunction

        function void connect_phase(uvm_phase phase);
            super.connect_phase(phase);
            drv.alu_vif = alu_cfg.alu_vif;
            mon.alu_vif = alu_cfg.alu_vif;
            drv.seq_item_port.connect(sqr.seq_item_export);
            mon.mon_ap.connect(agt_ap); // connect monitor share point with agent share point so the monitor will be able to get data from the scoreboard and the cov collector
        endfunction
    endclass
endpackage