import alu_seq_item_pkg::*;
class alu_coverage;

    mailbox #(alu_seq_item) mon2cov;

    alu_seq_item cov_item;

    // cover group
    covergroup alu_Cross_Group;
        // cover points
        A_CP: coverpoint cov_item.a {
            bins a_data_0 = {0};
            bins higher_half = {[8:14]};
            bins a_data_max = {15};
            bins a_data_default = default; 
        }

        B_CP: coverpoint cov_item.b {
            bins b_data_0 = {0};
            bins higher_half = {[8:14]};
            bins b_data_max = {15};
            bins b_data_default = default; 
        }

        OP_CP: coverpoint cov_item.op {
            bins op_add = {2'b00};
            bins op_xor  = {2'b01};
            bins op_and = {2'b10};
            bins op_or = {2'b11};
        }

        OUT_CP: coverpoint cov_item.out;
        C_CP: coverpoint cov_item.c;

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

    function new();
        alu_Cross_Group = new();
    endfunction

    function void connect_mail(mailbox #(alu_seq_item) m);
        mon2cov = m;
    endfunction

    task run();
        alu_seq_item item;
        forever begin
            mon2cov.get(item);
            cov_item = item;
            alu_Cross_Group.sample();
        end
    endtask

endclass


/*
package alu_coverage_pkg;
    import alu_seq_item_pkg::*;
    class alu_coverage;
        // alu_seq_item class object
        alu_seq_item alu_tr_el = new;

        // cover group
        covergroup alu_Cross_Group;
            // cover points
            A_CP: coverpoint a {
                bins a_data_0 = {0};
                bins higher_half = {[8:14]};
                bins a_data_max = {15};
                bins a_data_default = default; 
            }

            B_CP: coverpoint b {
                bins b_data_0 = {0};
                bins higher_half = {[8:14]};
                bins b_data_max = {15};
                bins b_data_default = default; 
            }

            OP_CP: coverpoint op {
                bins op_add = {2'b00};
                bins op_xor  = {2'b01};
                bins op_and = {2'b10};
                bins op_or = {2'b11};
            }

            OUT_CP: coverpoint out;
            C_CP: coverpoint c;

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

        // sample function
        function void sample_data(alu_seq_item alu_tr);
            alu_tr_el = alu_tr;
            alu_Cross_Group.sample();
        endfunction

        function new();
            alu_Cross_Group = new();
        endfunction
    endclass
endpackage*/