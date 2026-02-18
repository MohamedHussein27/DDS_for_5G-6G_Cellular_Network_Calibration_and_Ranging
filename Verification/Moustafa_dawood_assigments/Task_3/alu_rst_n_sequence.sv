package rst_sequence_pkg;

  import uvm_pkg::*;
  import alu_sequence_item_pkg::*;
  `include "uvm_macros.svh"

  class rst_sequence extends uvm_sequence #(alu_sequence_item);
    `uvm_object_utils(rst_sequence)

    function new(string name = "rst_sequence");
      super.new(name);
    endfunction


    task body();

      alu_sequence_item item;

      // ---------------------------------
      // 1️ ASSERT RESET
      // ---------------------------------
      repeat (3) begin
        item = alu_sequence_item::type_id::create("item");

        start_item(item);

        if (!item.randomize() with {
              rst_n == 0;
              a == 0;
              b == 0;
              op == 0;
            })
          `uvm_fatal("RST_SEQ", "Randomization failed during reset assert")

        finish_item(item);
      end

      item = alu_sequence_item::type_id::create("item");

      start_item(item);

      if (!item.randomize() with {
            rst_n == 1;
            a == 0;
            b == 0;
            op == 0;
          })
        `uvm_fatal("RST_SEQ", "Randomization failed during reset deassert")

      finish_item(item);

    endtask

  endclass

endpackage
