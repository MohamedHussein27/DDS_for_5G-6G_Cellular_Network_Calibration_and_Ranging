package random_op_pkg;
    import uvm_pkg::*;
	import alu_sequence_item_pkg::*;
	`include "uvm_macros.svh"

    class random_op extends uvm_sequence #(alu_sequence_item);
        `uvm_object_utils(random_op)
        function new(string name = "random_op");
            super.new(name);

        endfunction //new()


    task body();
    repeat (20000) 
    begin
        
            alu_sequence_item item;
            item = alu_sequence_item::type_id::create("item");
            start_item(item);
            assert(item.randomize());
            
            
            finish_item(item);
            end// repeat
    endtask
    endclass //random_op extends superClass
    
endpackage