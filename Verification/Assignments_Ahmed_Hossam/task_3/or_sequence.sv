`include "seq_item.sv"
class or_sequence extends uvm_sequence #(seq_item);
`uvm_object_utils(or_sequence)
  function new(string name = "or_sequence");
    super.new(name);
  endfunction

  task body();
    seq_item item;
    repeat (250)begin
    item = seq_item::type_id::create("item");
    start_item(item);
    assert (item.randomize()) else $fatal("Randomization failed");
    item.op = alu_pack::OR_op; // Set operation to OR
    finish_item(item);
    end
  endtask
endclass