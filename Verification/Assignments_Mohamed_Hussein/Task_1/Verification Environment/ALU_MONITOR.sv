import alu_seq_item_pkg::*;
class alu_monitor;

    virtual alu_if.MONITOR vif;
    mailbox #(alu_seq_item) mon2sb;
    mailbox #(alu_seq_item) mon2cov;

    function new();
    endfunction

    function void set_vif(virtual alu_if.MONITOR vif);
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
            @(posedge vif.clk);
            #1; // after posedge

            item = new();
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



/*package alu_monitor_pkg;
    import alu_seq_item_pkg::*;
    import alu_coverage_pkg::*;
    import alu_scoreboard_pkg::*;
    import shared_pkg::*;
    class alu_monitor;
        
        // virtual interface
        virtual alu_if.MONITOR vif;

        function new(virtual alu_if.MONITOR vif);
            this.vif = vif;
        endfunction

        //alu_seq_item item_mon = new; // seq_item object
        alu_coverage cov_mon = new; // coverage object
        alu_scoreboard scb_mon = new; // scoreboard object

        task sample(alu_seq_item item);
            #1;
            item.a   = vif.a;
            item.b   = vif.b;
            item.op  = vif.op;
            item.out = vif.out;
            item.c   = vif.c;
            /*$display(
                "[AFTER ][MONITOR] t=%0t : a=%0d b=%0d op=%0d | out=%0d c=%0d",
                $time, item.a, item.b, item.op, item.out, item.c
            );
            // two parallel processes
            fork
                // process 1
                cov_mon.sample_data(item);
                // process 2
                scb_mon.check_data(item);
            join
        endtask
    endclass
endpackage*/
        