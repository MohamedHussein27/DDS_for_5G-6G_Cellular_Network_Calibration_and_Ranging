package alu_driver_pkg;

import uvm_pkg::*;
import alu_pkg::*;
import alu_config_pkg::*;
	import alu_sequence_item_pkg::*;
	`include "uvm_macros.svh"

    class alu_driver extends uvm_driver #(alu_sequence_item);
        `uvm_component_utils(alu_driver)

        virtual alu_if alu_vif; // creating the virtual interface handle to drive the signls
        alu_config cfg;

        function new(string name = "alu_driver", uvm_component parent = null);
            super.new(name, parent);
        endfunction : new

    // need to revist this build phase
    
        function void build_phase(uvm_phase phase);
        super.build_phase(phase);
        //if(!uvm_config_db#(alu_config)::get(this,"","CFG",cfg))
       // `uvm_fatal("DRV","Cannot get config")
        // alu_vif = cfg.alu_vif;
        endfunction

        
        task run_phase(uvm_phase phase);
        alu_sequence_item item;
            super.run_phase(phase);
            
            forever begin
                item=alu_sequence_item::type_id::create("item");
                seq_item_port.get_next_item(item);
                if (item != null) begin
                    // drive the interface signals with the item data
                    alu_vif.op_set <= item.op;
                    alu_vif.operand_a <= item.a;
                    alu_vif.operand_b <= item.b;
                    alu_vif.rst_n <= item.rst_n;
                    // wait for one clock cycle
                    @(negedge alu_vif.clk);
                    // complete the item
                    seq_item_port.item_done();
                    `uvm_info("run_phase", item.convert2string(), UVM_HIGH)
                end
            end
        endtask : run_phase

    endclass




    
endpackage