import counter_seq_item_pkg::*;
class counter_monitor;
    string name;

    virtual counter_if counter_vif;
    mailbox #(counter_seq_item) mon2sb;
    mailbox #(counter_seq_item) mon2cov;

    function new(string name = "MONITOR");
        this.name = name;
        mon2sb = new();
        mon2cov = new();
    endfunction

    task run();
        counter_seq_item item;
        forever begin
            item = new();
            @(counter_vif.mon_cb); // wait for monitor clocking block to sample data
           
            item.rst_n          = counter_vif.rst_n;
            item.start          = counter_vif.start;
            item.wait_timer     = counter_vif.wait_timer;
            item.flag           = counter_vif.flag;
            item.busy           = counter_vif.busy;
            item.count_value    = counter_vif.count_value;

            // sending mailbox to scoreboard and coverage
            mon2sb.put(item);
            mon2cov.put(item);
            $display(
                "[AFTER ][MONITOR] t=%0t : start=%0d wait_timer=%0d flag=%0d | busy=%0d count_value=%0d",
                $time, item.start, item.wait_timer, item.flag, item.busy, item.count_value
            );
        end
    endtask

endclass
