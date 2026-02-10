// this module is to act as the sequence and sequencer in comparison with uvm
import counter_seq_item_pkg::*;
class counter_generator;
    string name;
    mailbox #(counter_seq_item) gen2drv;

    event gen_ack; // acknowledgment event to send transactions to driver

    bit has_started; // flag to indicate if start signal is rising edge or not
    bit has_reset; // flag to indicate if reset has been deasserted, so we will make start signal more frequent after reset is deasserted

    function new(string name = "counter_generator");
        this.name = name;
        gen2drv = new();
    endfunction

    task run();
        counter_seq_item item;
        item = new();
        // reset first (acting like reset sequence)
        item.rst_n      = 0;
        item.start      = 0;
        item.wait_timer = 0;
        item.flag       = 0;
        has_started     = item.start; // initialize has_started with the first transaction's start value
        has_reset       = item.rst_n; // initialize has_reset with the first transaction's rst_n value
        gen2drv.put(item);
        @(gen_ack);

        // directed tests (assigning the wait_timer value first)
        item = new();
        item.rst_n      = 1;
        item.start      = 0; 
        item.wait_timer = 16'h0003;
        has_started     = item.start;
        has_reset       = item.rst_n;
        gen2drv.put(item);
        @(gen_ack);


        item = new();
        item.rst_n      = 1;
        item.start      = 1; // start signal rising edge
        item.wait_timer = 16'h0003;
        has_started     = item.start;
        has_reset       = item.rst_n;
        gen2drv.put(item);
        @(gen_ack);

        // test 1 : testing wait_timer without flag
        for (int i = 0; i < 10; i++) begin
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = 16'h0003;
            has_started     = item.start;
            has_reset       = item.rst_n;
            gen2drv.put(item);
            @(gen_ack);
        end

        // wait for two transactions
        item = new();
        item.rst_n      = 1;
        item.wait_timer = 16'h0003;
        gen2drv.put(item);
        @(gen_ack);

        item = new();
        item.rst_n      = 1;
        item.wait_timer = 16'h0003;
        gen2drv.put(item);
        @(gen_ack);

        // test 2 : testing flag at internal counter != wait_timer
        for (int i = 0; i < 10; i++) begin
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = 16'h0003;
            if (i == 4) 
            item.flag       = 1;
            has_started     = item.start;
            has_reset       = item.rst_n;
            gen2drv.put(item);
            @(gen_ack);
        end

        item = new();
        item.rst_n      = 1;
        item.wait_timer = 16'h0003;
        gen2drv.put(item);
        @(gen_ack);

        // test 3 : testing flag at internal counter = wait_timer
        for (int i = 0; i < 10; i++) begin
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = 16'h0003;
            if (i == 4) 
            item.flag       = 1;
            has_started     = item.start;
            has_reset       = item.rst_n;
            gen2drv.put(item);
            @(gen_ack);
        end

        item = new();
        item.rst_n      = 1;
        item.wait_timer = 16'h0003;
        gen2drv.put(item);
        @(gen_ack);

        // test 4: testing the start without previously resetting the design
        item = new();
        item.rst_n      = 1;
        item.start      = 1;
        item.wait_timer = 16'h0003;
        gen2drv.put(item);
        @(gen_ack);

        for (int i = 0; i < 10; i++) begin
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = 16'h0003;
            item.flag       = 0;
            has_started     = item.start;
            has_reset       = item.rst_n;
            gen2drv.put(item);
            @(gen_ack);
        end

        // test 5: asserting start while normal operation
        item = new();
        item.rst_n      = 1;
        item.start      = 1;
        item.wait_timer = 16'h0003;
        gen2drv.put(item);
        @(gen_ack);

        for (int i = 0; i < 10; i++) begin
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = 16'h0003;
            item.flag       = 0;
            has_started     = item.start;
            has_reset       = item.rst_n;
            gen2drv.put(item);
            @(gen_ack);
        end

        //Test 6: testing max value of count
        item = new();
        item.rst_n      = 0; // reset first
        item.start      = 0; 
        item.wait_timer = 16'h0002;
        has_started     = item.start;
        has_reset       = item.rst_n;
        gen2drv.put(item);
        @(gen_ack);
        
        item = new();
        item.rst_n      = 1;
        item.start      = 1;
        item.wait_timer = 16'h0002;
        gen2drv.put(item);
        @(gen_ack);

        for (int i = 0; i < 70; i++) begin
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = 16'h0002;
            item.flag       = 0;
            has_started     = item.start;
            has_reset       = item.rst_n;
            gen2drv.put(item);
            @(gen_ack);
        end

        item = new();
        item.rst_n      = 1;
        item.start      = 0;
        item.wait_timer = 16'h0002;
        gen2drv.put(item);
        @(gen_ack);

        // test7 : testing asserting start after counter reaches max
        item = new();
        item.rst_n      = 1;
        item.start      = 1;
        item.wait_timer = 16'h0002;
        gen2drv.put(item);
        @(gen_ack);

        for (int i = 0; i < 10; i++) begin
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = 16'h0002;
            item.flag       = 0;
            has_started     = item.start;
            has_reset       = item.rst_n;
            gen2drv.put(item);
            @(gen_ack);
        end

        repeat (1000) begin
            item = new(); // every loop we create new sequence item (this is exactly like uvm_sequence_item::type_id::create() in sequence)
            if (!has_reset && !has_started) begin
                assert(item.randomize() with{ start dist {0 := 35, 1 := 65}; wait_timer == 2;}); // if reset is not deasserted, then it should be 1 in the next transaction
            end
            
            else if (has_started) begin
                assert(item.randomize() with { start == 0; wait_timer == 2;}); // if start is already 1, then it should be 0 in the next transaction
                has_started = 0;
                has_reset   = 1; // update has_reset to indicate that reset has been deasserted
            end
            else
                assert(item.randomize() with { wait_timer == 2; });
            
            has_started     = item.start; // update has_started with the current transaction's start value
            if (item.rst_n == 0) 
            has_reset       = 0; 
            
            gen2drv.put(item); // send the seq_item to driver via mailbox
            @(gen_ack); // wait for acknowledgment from driver before sending next transaction
        end
    endtask

endclass
