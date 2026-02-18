package or_sequence_pkg;

	import uvm_pkg::*;
	import alu_sequence_item_pkg::*;
	`include "uvm_macros.svh"

    class or_sequence extends uvm_sequence #(alu_sequence_item);
        `uvm_object_utils(or_sequence)

        function new(string name = "or_sequence");
            super.new(name);
        endfunction //new()


         task body();
            repeat (500) begin
            alu_sequence_item item;
            item = alu_sequence_item::type_id::create("item");
            start_item(item);
            assert(item.randomize() with { op == alu_pkg::OR_op; });
            // or we can use item.op = alu_pkg::OR_op; same result
            
            finish_item(item);
            end// repeat
        endtask //body()
    endclass
    
endpackage