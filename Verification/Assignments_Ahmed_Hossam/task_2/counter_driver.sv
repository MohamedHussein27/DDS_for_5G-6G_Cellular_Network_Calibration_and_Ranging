`include "counter_item.sv"
class counter_driver;
virtual counter_if vif;
mailbox #(counter_item) gen2driv;
     function new();
    endfunction //new()
    function void connecting (virtual counter_if vif , mailbox gen2driv);
    this.vif=vif;
    this.gen2driv=gen2driv;
    endfunction
    counter_item item;
    task run();
    forever begin
    item=new();
    gen2driv.get(item);
    vif.rst_n=item.rst_n;
    vif.start=item.start;
    vif.wait_timer=item.wait_timer;
    vif.flag=item.flag;
    @(negedge vif.clk);
    end
    endtask
endclass //counter_dri