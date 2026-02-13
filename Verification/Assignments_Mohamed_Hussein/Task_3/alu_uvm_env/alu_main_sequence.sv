package alu_main_sequence_pkg;
    import alu_seq_item_pkg::*;
    import alu_shared_pkg::*;
    import uvm_pkg::*;
    `include "uvm_macros.svh"
    class alu_main_sequence extends uvm_sequence #(alu_seq_item);
        `uvm_object_utils(alu_main_sequence);
        alu_seq_item seq_item;
        
        // constructor 
        function new(string name = "alu_main_sequence");
            super.new(name);
        endfunction

        task body;
            repeat (15000) begin
                seq_item = alu_seq_item::type_id::create("seq_item");
                start_item(seq_item);
                assert(seq_item.randomize());
                finish_item(seq_item);
                end
        endtask
    endclass
endpackage
