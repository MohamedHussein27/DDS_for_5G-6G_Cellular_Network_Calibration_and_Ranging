`include "seq_item.sv"
class driver extends uvm_driver #(seq_item);
`uvm_component_utils (driver)
    virtual alu_if vif;
    seq_item item;

    function new(string name = "driver", uvm_component parent = null);
        super.new(name, parent);
        endfunction
        
    task run_phase(uvm_phase phase);
    super.run_phase(phase);
        $display("[DRIVER] Waiting for transactions...");
        forever begin
            item =seq_item::type_id::create("item");
            seq_item_port.get_next_item(item); // Wait for item from Sequencer
            
           
            // Drive signals
            vif.a  = item.a;
            vif.b  = item.b;
            vif.op = item.op;
            @(negedge vif.clk);
            // Wait for propagation
            seq_item_port.item_done(); // Notify Sequencer that item is done
            `uvm_info("run_phase",item.convert2string(),UVM_MEDIUM)
            
        end
    endtask
endclass