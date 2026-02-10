`include "counter_generator.sv"
`include "counter_driver.sv"
`include "counter_monitor.sv"


class counter_agent;
counter_monitor mon;
counter_driver driv;
counter_generator gen;
virtual counter_if vif;
counter_coverage cov;  // The coverage component
mailbox #(counter_item) gen2driv;
mailbox #(counter_item) mon2scb;
    function new();
    endfunction //new()
    function void connecting (virtual counter_if vif);
    this.vif=vif;
    mon=new();
    driv=new();
    gen=new();
    cov=new();
    gen2driv=new();
    mon2scb=new();
    gen.connecting(gen2driv);
    driv.connecting(vif,gen2driv);
    mon.connecting(vif,mon2scb,cov);
        
    endfunction
    task run();
    fork
        driv.run();
        mon.run();
        
    join_none
    gen.run();
    endtask

endclass //counter_env