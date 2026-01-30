package alu_seq_item_pkg;
    import shared_pkg::*;
    class alu_seq_item;
        rand bit [3:0] a;
        rand bit [3:0] b;
        rand bit [1:0] op;
        logic [3:0] out;
        logic c;

        logic [7:0] temp_address; // will be used is equal address constraint

        // constraints
        constraint a_con {
            a dist {[0:12] := 70, [13:15] := 30}; // higher values has less prop.
        }

        constraint b_con {
            b dist {[0:14] := 70, 15 := 30}; // higher values has less prop.
        }
        
    endclass
endpackage
                
            