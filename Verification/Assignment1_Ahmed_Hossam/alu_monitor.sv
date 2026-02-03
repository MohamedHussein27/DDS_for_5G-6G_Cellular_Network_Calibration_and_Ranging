`include "alu_item.sv"
`include "alu_coverage.sv"
class alu_monitor;
    virtual alu_if vif;
    mailbox #(alu_item) mon2scb;
    alu_coverage cov;
    function new(virtual alu_if vif, mailbox mon2scb,alu_coverage cov);
        this.vif = vif;
        this.mon2scb = mon2scb;
        this.cov = cov; // Store the handle
    endfunction

    task run();
        forever begin
            alu_item item = new(); // Always create a NEW object
            
           @(posedge vif.clk);
           #1; // Small delay to avoid
            
            // Sample
            item.a = vif.a;
            item.b = vif.b;
            void'($cast(item.op, vif.op)); // Cast 2-bit wire to enum
            item.c = vif.c;
            item.out = vif.out;

            cov.sample(item);
            mon2scb.put(item);

        end
    endtask
endclass