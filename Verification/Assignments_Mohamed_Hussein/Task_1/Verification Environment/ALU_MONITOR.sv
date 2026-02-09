import alu_seq_item_pkg::*;
class alu_monitor;

    virtual alu_if vif;
    mailbox #(alu_seq_item) mon2sb;
    mailbox #(alu_seq_item) mon2cov;

    event mon_start; // monitor start event to synchronize with driver

    function new();
    endfunction

    function void set_vif(virtual alu_if vif);
        this.vif = vif;
    endfunction

    function void connect_mails(
            mailbox #(alu_seq_item) mon2sb,
            mailbox #(alu_seq_item) mon2cov
        );
        this.mon2sb  = mon2sb;
        this.mon2cov = mon2cov;
    endfunction

    task run();
        alu_seq_item item;
        forever begin
            item = new();
            @(posedge vif.clk);
            @(mon_start); // synchronize with driver

            
            item.rst_n = vif.rst_n;
            item.a     = vif.a;
            item.b     = vif.b;
            item.op    = vif.op;
            item.out   = vif.out;
            item.c     = vif.c;

            // sending mailbox to scoreboard and coverage
            mon2sb.put(item);
            mon2cov.put(item);
            $display(
                "[AFTER ][MONITOR] t=%0t : a=%0d b=%0d op=%0d | out=%0d c=%0d",
                $time, item.a, item.b, item.op, item.out, item.c
            );
        end
    endtask

endclass

        