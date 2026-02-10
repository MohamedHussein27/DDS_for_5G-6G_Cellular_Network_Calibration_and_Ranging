`include "counter_item.sv"
`include "counter_coverage.sv"
class counter_monitor;
virtual counter_if vif;
mailbox #(counter_item) mon2scb;
counter_coverage cov;
    function new();
    endfunction //new()
    function void connecting (virtual counter_if vif , mailbox mon2scb,counter_coverage cov);
    this.vif=vif;
    this.mon2scb=mon2scb;
    this.cov=cov;
    endfunction
    task run ();
    forever begin
        counter_item item;
        item=new();
        @(posedge vif.clk);
        #1 ;
        item.rst_n=vif.rst_n;
        item.start=vif.start;
        item.wait_timer=vif.wait_timer;
        item.flag=vif.flag;
        item.count_value=vif.count_value;
        item.busy=vif.busy;
        cov.sample(item);
        mon2scb.put(item);

    
    end
        
    endtask //
endclass //counter_mo