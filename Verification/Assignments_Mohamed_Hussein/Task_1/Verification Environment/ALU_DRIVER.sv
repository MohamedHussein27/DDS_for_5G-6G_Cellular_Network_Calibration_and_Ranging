import alu_seq_item_pkg::*;
class alu_driver;

    virtual alu_if vif;
    mailbox #(alu_seq_item) gen2drv;

    event drv_rqt; // driver request event to acknowledge generator

    function new();
    endfunction

    function void set_vif(virtual alu_if vif);
        this.vif = vif;
    endfunction

    function void connect_mail(mailbox #(alu_seq_item) m);
        gen2drv = m;
    endfunction

    task run();
        alu_seq_item item;
        forever begin
            item = new();
            gen2drv.get(item);
            @(negedge vif.clk);

            vif.rst_n = item.rst_n;
            vif.a     = item.a;
            vif.b     = item.b;
            vif.op    = item.op;
            $display(
                "[BEFORE][DRIVER] t=%0t : a=%0d b=%0d op=%0d",
                $time, item.a, item.b, item.op
            );
            
            -> drv_rqt; // notify generator that transaction is done
        end
    endtask

endclass
