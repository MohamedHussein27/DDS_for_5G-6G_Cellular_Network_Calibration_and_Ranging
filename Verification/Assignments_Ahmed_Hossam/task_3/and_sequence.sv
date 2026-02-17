`include "seq_item.sv"
class and_sequence extends uvm_sequence #(seq_item);
`uvm_object_utils(and_sequence)
  function new(string name = "and_sequence");
    super.new(name);
  endfunction

  task body();
    seq_item item;
    repeat (250)begin
    item = seq_item::type_id::create("item");
    start_item(item);
    assert (item.randomize()) else $fatal("Randomization failed");
    item.op = alu_pack::AND_op; // Set operation to AND
    finish_item(item);
    end
  endtask
endclass