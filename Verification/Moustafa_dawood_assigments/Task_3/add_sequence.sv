package add_sequence_pkg;

	import uvm_pkg::*;
	import alu_sequence_item_pkg::*;
	`include "uvm_macros.svh"

    class add_sequence extends uvm_sequence #(alu_sequence_item);
        `uvm_object_utils(add_sequence)

        function new(string name = "add_sequence");
            super.new(name);
        endfunction //new()

// virtual task or not ??
         task body();
            repeat (500) begin
            alu_sequence_item item;
            item = alu_sequence_item::type_id::create("item");
            start_item(item);
            assert(item.randomize() with { op == alu_pkg::ADD_op; });
            // or we can use item.op = alu_pkg::ADD_op; same result
            
            finish_item(item);
            end// repeat
        endtask //body()
    endclass
    
endpackage