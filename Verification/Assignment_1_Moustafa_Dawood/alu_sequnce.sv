/*
class alu_sequence;

  mailbox #(alu_seq_item) seq2drv_mb;

  function new(mailbox #(alu_seq_item) mb);
    seq2drv_mb = mb;
  endfunction

  task run();
    alu_seq_item item;

    repeat (10000) begin
      item = new();
      assert(item.randomize());
      seq2drv_mb.put(item);
    end
  endtask

endclass
*/


class alu_sequence;

  mailbox #(alu_seq_item) seq2drv_mb;

  function new(mailbox #(alu_seq_item) mb);
    seq2drv_mb = mb;
  endfunction

  task run();
    alu_seq_item item;

    for (int op = 0; op < 4; op++) begin
      for (int a = 0; a < 16; a++) begin
        for (int b = 0; b < 16; b++) begin
          item = new();
          item.op = opcode_t'(op);  //  casting int to enum
          item.a  = a;
          item.b  = b;
          seq2drv_mb.put(item);
        end
      end
    end
  endtask

endclass


