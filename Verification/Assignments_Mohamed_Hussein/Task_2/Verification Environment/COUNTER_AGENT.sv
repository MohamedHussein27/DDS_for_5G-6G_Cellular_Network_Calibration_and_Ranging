`include "COUNTER_GENERATOR.sv"
`include "COUNTER_DRIVER.sv"
`include "COUNTER_MONITOR.sv"
import counter_seq_item_pkg::*;

class counter_agent;
    string name;

    counter_generator gen;
    counter_driver    drv;
    counter_monitor   mon;

    // virtual interface
    virtual counter_if counter_vif;


    mailbox #(counter_seq_item) mon2sb;
    mailbox #(counter_seq_item) mon2cov;

    function new(string name = "counter_agent");
        this.name = name;
        mon2sb = new();
        mon2cov = new();
        gen = new();
        drv = new();
        mon = new();
    endfunction

    function void connect();
        // connect agent virtual interface to driver and monitor virtual interface
        drv.counter_vif = counter_vif;
        mon.counter_vif = counter_vif;
        // connect mailboxs to monitor, generator and driver
        mon.mon2sb  = mon2sb;
        mon.mon2cov =  mon2cov;
        drv.gen2drv = gen.gen2drv;
        gen.gen_ack = drv.drv_rqt; // connect acknowledgment event (TLM (port, export) in UVM)
    endfunction

    task run();
        fork
            gen.run();
            drv.run();
            mon.run();
        join_any
    endtask

endclass
