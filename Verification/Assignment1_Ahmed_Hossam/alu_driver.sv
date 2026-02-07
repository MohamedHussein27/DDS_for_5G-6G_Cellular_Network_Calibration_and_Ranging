`include "alu_item.sv"
class alu_driver;
    virtual alu_if vif;
    mailbox #(alu_item) gen2driv;
    function new();
        endfunction
        
    function void connecting (virtual alu_if vif, mailbox gen2driv);
        this.vif = vif;
        this.gen2driv = gen2driv;
    endfunction

    task run();
        $display("[DRIVER] Waiting for transactions...");
        forever begin
            alu_item item;
            @(negedge vif.clk);
            gen2driv.get(item); // Wait for item from Generator
            
            // Drive signals
            vif.a  = item.a;
            vif.b  = item.b;
            vif.op = item.op;
            
            // Wait for propagation
            #10; 
        end
    endtask
endclass