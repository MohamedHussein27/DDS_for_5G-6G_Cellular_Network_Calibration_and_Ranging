`include "ALU_AGENT.sv"
`include "ALU_COVERAGE.sv"
`include "ALU_SCOREBOARD.sv"
import shared_pkg::*;
class alu_env;

    // virtual interface
    virtual alu_if alu_vif;

    // agent, scoreboard and coverage class objects
    alu_agent agt;
    alu_scoreboard sb;
    alu_coverage cov;

    // mailbox for scoreboard and coverage
    mailbox #(alu_seq_item) mon2sb;
    mailbox #(alu_seq_item) mon2cov;

    function new();
        agt = new();
        sb  = new();
        cov = new();
        mon2sb = new();
        mon2cov = new();
    endfunction

    function void set_vif(virtual alu_if alu_vif);
        this.alu_vif = alu_vif;
        agt.set_vif_drv(alu_vif.TEST);
        agt.set_vif_mon(alu_vif.MONITOR);
    endfunction

    // function to connect seq_item to scoreboard and coverage
    function void connect_mails();
        agt.connect_mails(mon2sb, mon2cov);
        sb.connect_mail(mon2sb);
        cov.connect_mail(mon2cov);
    endfunction

    task run();
        fork
            agt.run();
            sb.run();
            cov.run();
        join_none
    endtask

endclass