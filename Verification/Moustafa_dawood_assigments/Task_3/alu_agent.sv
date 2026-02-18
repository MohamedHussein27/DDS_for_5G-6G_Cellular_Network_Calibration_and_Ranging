package alu_agent_pkg;
	import uvm_pkg::*;
    import alu_sequencer_pkg::*;
    import alu_driver_pkg::*;
    import alu_monitor_pkg::*;
    import alu_config_pkg::*;
    import alu_sequence_item_pkg::*;

    `include "uvm_macros.svh"

    class alu_agent extends uvm_agent;
        `uvm_component_utils(alu_agent)

        alu_sequencer sqr;
        alu_driver driver;
        alu_monitor monitor;
        alu_config alu_cfg;
        uvm_analysis_port#(alu_sequence_item) agt_ap; // analysis port for sending

        function new(string name = "alu_agent", uvm_component parent = null);
            super.new(name, parent);
        endfunction : new

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);
            
            if(!uvm_config_db#(alu_config)::get(this, "", "CFG", alu_cfg)) //get the cfg from the db and assign it to the cfg of the driver
				`uvm_fatal("build_phase", "Agent - Unable to get the configuration object")

            sqr = alu_sequencer::type_id::create("sqr", this);
            driver = alu_driver::type_id::create("driver", this);
            monitor = alu_monitor::type_id::create("monitor", this);

            agt_ap = new("agt_ap", this);

        endfunction : build_phase

        function void connect_phase(uvm_phase phase);

            super.connect_phase(phase);
            
         driver.alu_vif = alu_cfg.alu_vif; // the agent connect the vif of the driver to the vif of the cfg
          monitor.alu_vif = alu_cfg.alu_vif;

            driver.seq_item_port.connect(sqr.seq_item_export); // buit in connection betwenn the sqr and driver no need to create any ports
            //monitor.mon_ap.connect(sequencer.mon_port);
            monitor.mon_ap.connect(agt_ap);



        endfunction : connect_phase

    endclass
endpackage