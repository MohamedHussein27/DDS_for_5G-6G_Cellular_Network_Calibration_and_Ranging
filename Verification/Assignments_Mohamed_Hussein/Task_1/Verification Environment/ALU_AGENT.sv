`include "ALU_GENERATOR.sv"
`include "ALU_DRIVER.sv"
`include "ALU_MONITOR.sv"
import alu_seq_item_pkg::*;

class alu_agent;

    alu_generator gen;
    alu_driver    drv;
    alu_monitor   mon;

    mailbox #(alu_seq_item) gen2drv;
    mailbox #(alu_seq_item) mon2sb;
    mailbox #(alu_seq_item) mon2cov;

    function new();
        gen2drv = new();
        gen = new();
        drv = new();
        mon = new();
    endfunction

    function void set_vif_drv(virtual alu_if.TEST vif);
        drv.set_vif(vif);
    endfunction

    function void set_vif_mon(virtual alu_if.MONITOR vif);
        mon.set_vif(vif);
    endfunction

    function void connect_mails(
        mailbox #(alu_seq_item) mon2sb,
        mailbox #(alu_seq_item) mon2cov
    );
        this.mon2sb  = mon2sb;
        this.mon2cov = mon2cov;
        mon.connect_mails(mon2sb, mon2cov);
        drv.connect_mail(gen2drv);
        gen.connect_mail(gen2drv);
    endfunction

    task run();
        fork
            gen.run();
            drv.run();
            mon.run();
        join_none
    endtask

endclass
