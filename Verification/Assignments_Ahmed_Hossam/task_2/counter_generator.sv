`include "counter_item.sv"
class counter_generator ;
mailbox #(counter_item) gen2driv;

    function new();
    endfunction //new()
    function void connecting (mailbox gen2driv);
    this.gen2driv=gen2driv;
    endfunction
task run();
counter_item item;
item=new();
item.rst_n=0;
item.start=0;
item.flag=0;
gen2driv.put(item);
repeat (1000) begin
item =new();
assert (item.randomize());
gen2driv.put(item);
end
endtask
endclass //generator extends superClass