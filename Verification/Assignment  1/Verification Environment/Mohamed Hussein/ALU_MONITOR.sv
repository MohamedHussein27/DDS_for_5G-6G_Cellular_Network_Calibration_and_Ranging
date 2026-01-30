import alu_seq_item_pkg::*;
import alu_coverage_pkg::*;
import alu_scoreboard_pkg::*;
import shared_pkg::*;
module alu_monitor (alu_if.MONITOR aluif);
    alu_seq_item tr_mon = new; // seq_itemn object
    alu_coverage cov_mon = new; // coverage object
    alu_scoreboard scb_mon = new; // scoreboard object

    initial begin 
        forever begin
            #1; // negedge
            // assigning interface data to class seq_itemn object
            tr_mon.a = aluif.a;
            tr_mon.b = aluif.b;
            tr_mon.op = aluif.op;
            tr_mon.out = aluif.out;
            tr_mon.c = aluif.c;
            // two parallel processes
            fork
                // process 1
                cov_mon.sample_data(tr_mon);
                // process 2
                scb_mon.check_data(tr_mon);
            join
            if (test_finished) begin
                $display("error count = %0d, correct count = %0d", error_count_out, correct_count_out);
                $stop;
            end
        end
    end
endmodule
        