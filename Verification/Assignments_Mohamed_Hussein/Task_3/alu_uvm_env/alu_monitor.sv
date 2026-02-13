package alu_monitor_pkg;
    import alu_seq_item_pkg::*;
    import uvm_pkg::*;
    `include "uvm_macros.svh"   
    class alu_monitor extends uvm_monitor;
        `uvm_component_utils(alu_monitor)

        virtual alu_if alu_vif; // virtual interface
        alu_seq_item rsp_seq_item; // sequence item
        uvm_analysis_port #(alu_seq_item) mon_ap;

        function new(string name = "alu_monitor", uvm_component parent = null);
            super.new(name, parent);
        endfunction

        // building share point
        function void build_phase(uvm_phase phase);
            super.build_phase(phase);
            mon_ap = new("mon_ap", this);
        endfunction

        task run_phase(uvm_phase phase);
            super.run_phase(phase);
            forever begin
                rsp_seq_item = alu_seq_item::type_id::create("rsp_seq_item");
                @(negedge alu_vif.clk);
                
                rsp_seq_item.rst_n  = alu_vif.rst_n;
                rsp_seq_item.a      = alu_vif.a;
                rsp_seq_item.b      = alu_vif.b;
                rsp_seq_item.op     = alu_vif.op;
                rsp_seq_item.out    = alu_vif.out;
                rsp_seq_item.c      = alu_vif.c;

                mon_ap.write(rsp_seq_item);
                `uvm_info("run_phase", rsp_seq_item.convert2string(), UVM_HIGH);
            end
        endtask
    endclass
endpackage