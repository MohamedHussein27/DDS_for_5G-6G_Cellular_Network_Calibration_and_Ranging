package alu_xor_sequence_pkg;
    import alu_seq_item_pkg::*;
    import alu_shared_pkg::*;
    import uvm_pkg::*;
    `include "uvm_macros.svh"
    class alu_xor_sequence extends uvm_sequence #(alu_seq_item);
        `uvm_object_utils(alu_xor_sequence);
        alu_seq_item seq_item;
        
        // constructor 
        function new(string name = "alu_xor_sequence");
            super.new(name);
        endfunction

        task body;
            seq_item = alu_seq_item::type_id::create("seq_item");
            start_item(seq_item);
            seq_item.rst_n = 1;
            seq_item.a     = $random();
            seq_item.b     = $random();
            seq_item.op    = 1;
            finish_item(seq_item);
        endtask
    endclass
endpackage
