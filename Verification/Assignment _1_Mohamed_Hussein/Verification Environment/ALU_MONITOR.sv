package alu_monitor_pkg;
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

        //alu_seq_item tr_mon = new; // seq_item object
        alu_coverage cov_mon = new; // coverage object
        alu_scoreboard scb_mon = new; // scoreboard object

        task sample(ref alu_seq_item tr);
            #1;
            tr.a   = vif.a;
            tr.b   = vif.b;
            tr.op  = vif.op;
            tr.out = vif.out;
            tr.c   = vif.c;
            /*$display(
                "[AFTER ][MONITOR] t=%0t : a=%0d b=%0d op=%0d | out=%0d c=%0d",
                $time, tr.a, tr.b, tr.op, tr.out, tr.c
            );*/
            // two parallel processes
            fork
                // process 1
                cov_mon.sample_data(tr);
                // process 2
                scb_mon.check_data(tr);
            join
        endtask
    endclass
endpackage
        