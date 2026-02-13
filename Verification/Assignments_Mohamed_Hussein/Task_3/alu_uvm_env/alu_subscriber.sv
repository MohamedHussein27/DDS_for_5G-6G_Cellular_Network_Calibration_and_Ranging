package alu_subscriber_pkg;
    import alu_seq_item_pkg::*;
    import uvm_pkg::*;
    `include "uvm_macros.svh"
    class alu_subscriber extends uvm_component;
        `uvm_component_utils(alu_subscriber)
        uvm_analysis_export #(alu_seq_item) sub_export;
        uvm_tlm_analysis_fifo #(alu_seq_item) sub_fifo;

        alu_seq_item alu_item_cov;
        virtual alu_internal_if alu_internal_vif; // virtual interface (internal) (to monitor internals)

        // cover group
        covergroup alu_Group;
            // cover points
            A_CP: coverpoint alu_item_cov.a {
                bins a_data_0 = {0};
                bins higher_half = {[8:14]};
                bins a_data_max = {15};
                bins a_data_default = default; 
            }

            B_CP: coverpoint alu_item_cov.b {
                bins b_data_0 = {0};
                bins higher_half = {[8:14]};
                bins b_data_max = {15};
                bins b_data_default = default; 
            }

            OP_CP: coverpoint alu_item_cov.op {
                bins op_add = {2'b00};
                bins op_xor  = {2'b01};
                bins op_and = {2'b10};
                bins op_or = {2'b11};
            }

            OUT_CP: coverpoint alu_item_cov.out;
            C_CP: coverpoint alu_item_cov.c;

            // cross coverage
            BOUNDARY_ADD_C: cross A_CP, B_CP, OP_CP {
                bins ADD_BOUNDARY = binsof(A_CP) intersect {15} && 
                                    binsof(B_CP) intersect {15} && 
                                    binsof(OP_CP) intersect {2'b00};
                option.cross_auto_bin_max = 0;
            } 

            CARRY_C: cross A_CP, B_CP, OP_CP {
                bins CARRY_SIG = binsof(A_CP.higher_half) && 
                                    binsof(B_CP.higher_half) && 
                                    binsof(OP_CP) intersect {2'b00};
                option.cross_auto_bin_max = 0;
            }
        endgroup

        function new (string name = "alu_subscriber", uvm_component parent = null);
            super.new(name, parent);
            alu_Group = new;
        endfunction

        function void build_phase(uvm_phase phase);
            super.build_phase(phase);
            sub_export = new("sub_export", this);
            sub_fifo = new("sub_fifo", this);
        endfunction

        function void connect_phase(uvm_phase phase);
            super.connect_phase(phase);
            sub_export.connect(sub_fifo.analysis_export);
        endfunction

        task run_phase(uvm_phase phase);
            super.run_phase(phase);
            forever begin
                sub_fifo.get(alu_item_cov);        
                alu_Group.sample();
            end
        endtask
    endclass
endpackage
