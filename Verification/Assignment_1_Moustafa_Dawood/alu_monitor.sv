
class alu_monitor;

  mailbox #(alu_seq_item) drv2mon_mb;
  mailbox #(alu_seq_item) mon2sb_mb;

  function new(
    mailbox #(alu_seq_item) d2m,
    mailbox #(alu_seq_item) m2s
  );
    drv2mon_mb = d2m;
    mon2sb_mb  = m2s;
  endfunction

  task run();
    alu_seq_item item;
    forever begin
      drv2mon_mb.get(item);
      mon2sb_mb.put(item);
    end
  endtask

endclass
