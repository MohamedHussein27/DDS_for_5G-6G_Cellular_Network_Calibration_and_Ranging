import counter_seq_item_pkg::*;
import shared_pkg::*;
class counter_driver;
    string name;

    virtual counter_if counter_vif;
    mailbox #(counter_seq_item) gen2drv;

    event drv_rqt; // driver request event to acknowledge generator

    function new(string name = "DRIVER");
        this.name = name;
        gen2drv = new();
    endfunction

    task run();
        counter_seq_item item;
        forever begin
            item = new();
            gen2drv.get(item);
            @(negedge counter_vif.clk);
            counter_vif.flag        = item.flag;

            if (!stop_driving) begin
                
                counter_vif.rst_n       = item.rst_n;
                counter_vif.start       = item.start;
                counter_vif.wait_timer  = item.wait_timer;
                
                /*$display(
                    "[BEFORE][DRIVER] t=%0t : start=%0d wait_timer=%0d flag=%0d",
                    $time, item.start, item.wait_timer, item.flag
                );*/
            end
            -> drv_rqt; // notify generator that transaction is done
        end
    endtask

endclass
