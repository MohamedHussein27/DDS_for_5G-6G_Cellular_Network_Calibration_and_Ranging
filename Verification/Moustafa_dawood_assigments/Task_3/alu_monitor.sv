package alu_monitor_pkg;
    import uvm_pkg::*;
    import alu_pkg::*;
    import alu_config_pkg::*;

	import alu_sequence_item_pkg::*;
	`include "uvm_macros.svh"

    class alu_monitor extends uvm_monitor;
    `uvm_component_utils(alu_monitor)
     alu_sequence_item item;
     virtual alu_if alu_vif;
     alu_config alu_cfg;

    uvm_analysis_port #(alu_sequence_item) mon_ap; // creating the analysis prort to  broadcast data to subscribers 

      function new(string name = "alu_monitor", uvm_component parent = null);
			super.new(name, parent);
		endfunction 

    


    function void build_phase(uvm_phase phase);
			super.build_phase(phase);
			mon_ap = new("mon_ap", this);
           // if(!uvm_config_db#(alu_config)::get(this, "", "CFG", alu_cfg))
   //`uvm_fatal("MON", "Cannot get config")

   // alu_vif = alu_cfg.alu_vif;

		endfunction : build_phase



    task run_phase(uvm_phase phase);
        super.run_phase(phase);
        forever begin
            @(posedge alu_vif.clk); // wait for the posative edge of the clock to sample the signals

            #1; // small delay to ensure stable values
            //(we could use @negedge clk withould using posedge and delay)
            
            item = alu_sequence_item::type_id::create("item");
            item.op = alu_vif.op_set;
            item.a = alu_vif.operand_a;
            item.b = alu_vif.operand_b;
            item.rst_n = alu_vif.rst_n;
            item.result = alu_vif.result;
            `uvm_info("alu_monitor", item.convert2string(), UVM_HIGH)
            mon_ap.write(item); // broadcast the item to subscribers
            `uvm_info("run_phase", item.convert2string(), UVM_HIGH)
        end
    endtask

    endclass //className extends superClass

    
endpackage