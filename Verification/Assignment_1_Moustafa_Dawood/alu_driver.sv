

class alu_driver;
  

  virtual alu_bfm bfm;
  mailbox #(alu_seq_item) seq2drv_mb;
  mailbox #(alu_seq_item) drv2mon_mb;

  function new(
    virtual alu_bfm b,
    mailbox #(alu_seq_item) s2d,
    mailbox #(alu_seq_item) d2m
  );
    bfm = b;
    seq2drv_mb = s2d;
    drv2mon_mb = d2m;
  endfunction

  task run();
    alu_seq_item item;

    forever begin
      seq2drv_mb.get(item);

      bfm.operand_a = item.a;
      bfm.operand_b = item.b;
      bfm.op_set    = item.op;

      #1; // combinational settle
      item.result = bfm.result;

      drv2mon_mb.put(item);
    end
  endtask

endclass
