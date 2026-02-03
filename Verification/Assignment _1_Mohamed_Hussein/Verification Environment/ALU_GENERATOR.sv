// this module is to act as the sequence and sequencer in comparison with uvm
import alu_seq_item_pkg::*;
class alu_generator;

    mailbox #(alu_seq_item) gen2drv;

    function new();
    endfunction

    function void connect_mail(mailbox #(alu_seq_item) m);
        gen2drv = m;
    endfunction

    task run();
        alu_seq_item item;
        item = new();
        // reset first (acting like reset sequence)
        item.rst_n = 0;
        item.a     = 0;
        item.b     = 0;
        item.op    = 0;
        gen2drv.put(item);
        item = new();
        item.rst_n = 1;
        gen2drv.put(item);
        repeat (15000) begin
            item = new(); // every loop we create new sequence item (this is exactly like uvm_sequence_item::type_id::create() in sequence)
            assert(item.randomize());
            gen2drv.put(item); // send the seq_item to driver via mailbox
        end
    endtask

endclass
