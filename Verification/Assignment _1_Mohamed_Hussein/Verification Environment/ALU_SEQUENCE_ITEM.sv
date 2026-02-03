
import shared_pkg::*;
package alu_seq_item_pkg;
    class alu_seq_item;
        bit clk;
        rand bit rst_n;
        rand bit [3:0] a;
        rand bit [3:0] b;
        rand bit [1:0] op;
        logic [3:0] out;
        logic c;

        // constraints
        constraint rst_n_con {
            rst_n dist {0 := 1, 1 := 99}; // reset active low
        }

        constraint a_con {
            a dist {[0:12] := 70, [13:15] := 30}; // higher values has less prop.
        }

        constraint b_con {
            b dist {[0:14] := 70, 15 := 30}; // higher values has less prop.
        }
        
    endclass
endpackage
            