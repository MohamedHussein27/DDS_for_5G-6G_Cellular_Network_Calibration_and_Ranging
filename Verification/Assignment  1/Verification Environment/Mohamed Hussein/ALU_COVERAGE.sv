package alu_coverage_pkg;
    import alu_seq_item_pkg::*;
    class alu_coverage;
        // alu_seq_item class object
        alu_seq_item alu_tr_el = new;

        // cover group
        covergroup alu_Cross_Group;
            // cover points
            A_CP: coverpoint alu_tr_el.a {
                bins a_data_0 = {0};
                bins higher_half = {[8:14]};
                bins a_data_max = {15};
                bins a_data_default = default; 
            }

            B_CP: coverpoint alu_tr_el.b {
                bins b_data_0 = {0};
                bins higher_half = {[8:14]};
                bins b_data_max = {15};
                bins b_data_default = default; 
            }

            OP_CP: coverpoint alu_tr_el.op {
                bins op_add = {2'b00};
                bins op_xor  = {2'b01};
                bins op_and = {2'b10};
                bins op_or = {2'b11};
            }

            OUT_CP: coverpoint alu_tr_el.out;
            C_CP: coverpoint alu_tr_el.c;

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
endpackage