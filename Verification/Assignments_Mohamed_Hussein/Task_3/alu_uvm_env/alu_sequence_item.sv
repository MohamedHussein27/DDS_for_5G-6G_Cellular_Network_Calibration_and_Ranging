package alu_seq_item_pkg;
    import uvm_pkg::*;
    `include "uvm_macros.svh"

    class alu_seq_item extends uvm_sequence_item;
        `uvm_object_utils(alu_seq_item)

        rand bit rst_n;
        rand bit [3:0] a;
        rand bit [3:0] b;
        rand bit [1:0] op;
        logic [3:0] out;
        logic c;

        function new(string name = "alu_seq_item");
            super.new(name);
        endfunction

        function string convert2string();
            return $sformatf("%s rst_n = 0b%0b, A = 0b%0b, B = 0b%0b, OP = 0b%0b, OUT = 0b%0b, Carry = 0b%0b",
            super.convert2string(), rst_n, a, b, op, out, c);
        endfunction

        function string convert2string_stimulus();
            return $sformatf("%s rst_n = 0b%0b, A = 0b%0b, B = 0b%0b, OP = 0b%0b",
            super.convert2string(), rst_n, a, b, op);
        endfunction

        constraint rst_n_con {
            rst_n dist {0 := 1, 1 := 99};
        }

        constraint a_con {
            a dist {[0:12] := 70, [13:15] := 30};
        }

        constraint b_con {
            b dist {[0:14] := 70, 15 := 30};
        }
    endclass
endpackage