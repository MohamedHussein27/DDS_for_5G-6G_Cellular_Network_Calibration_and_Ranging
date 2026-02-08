import alu_seq_item_pkg::*;
class alu_driver;

    virtual alu_if.TEST vif;
    mailbox #(alu_seq_item) gen2drv;

    function new();
    endfunction

    function void set_vif(virtual alu_if.TEST vif);
        this.vif = vif;
    endfunction

    function void connect_mail(mailbox #(alu_seq_item) m);
        gen2drv = m;
    endfunction

    task run();
        alu_seq_item item;
        forever begin
            gen2drv.get(item);

            vif.rst_n = item.rst_n;
            vif.a     = item.a;
            vif.b     = item.b;
            vif.op    = item.op;
            $display(
                "[BEFORE][DRIVER] t=%0t : a=%0d b=%0d op=%0d",
                $time, item.a, item.b, item.op
            );
            @(negedge vif.clk);
        end
    endtask

endclass






/*package alu_driver_pkg;
    import alu_seq_item_pkg::*;
    import shared_pkg::*;
    class alu_driver;

        virtual alu_if.TEST vif; // virtual interface
        
        function new(virtual alu_if.TEST vif);
            this.vif = vif;
        endfunction 

        // class object
        //alu_seq_item seq_dr = new;

        // stimulus
        task drive(alu_seq_item tr);
            vif.a  = tr.a;
            vif.b  = tr.b;
            vif.op = tr.op;
            /*$display(
                "[BEFORE][DRIVER] t=%0t : a=%0d b=%0d op=%0d",
                $time, tr.a, tr.b, tr.op
            );
            #1; // settle time for combinational dut
        endtask
    endclass
endpackage*/