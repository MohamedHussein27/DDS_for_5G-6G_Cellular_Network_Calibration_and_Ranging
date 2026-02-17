// this module is to act as the sequence and sequencer in comparison with uvm
import counter_seq_item_pkg::*;
import shared_pkg::*;
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
        
        // outer loop
        for (int j = 1; j < 301; j++) begin
            item = new();
            // reset first (acting like reset sequence)
            item.rst_n      = 0;
            item.start      = 0;
            item.wait_timer = 0;
            item.flag       = 0;            
            gen2drv.put(item);
            @(gen_ack);

            // directed tests (assigning the wait_timer value first)
            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.wait_timer = $random();
            waiting_time    = item.wait_timer;
            gen2drv.put(item);
            @(gen_ack);
            
            item = new();
            item.rst_n      = 1;
            item.start      = 1; // start signal rising edge
            item.wait_timer = waiting_time;
            waiting_time    = waiting_time;
            inner1 = 0;
            inner2 = 0;
            inner3 = 0;
            gen2drv.put(item);
            @(gen_ack);

            // testing normal operation with no flag
            for ( i = 0; i < (waiting_time*10 + 1); i++) begin
                item = new();
                gen2drv.put(item);
                @(gen_ack);
                inner1 = 1;
                inner2 = 0;
                inner3 = 0;
            end

            // testing flag while not wait_timer interval edge (glitch)
            if (waiting_time != 1) begin
                for ( i = 0; i < (waiting_time*10); i++) begin
                    item = new();
                    if (i >= waiting_time*5 && waiting_time == 2) begin
                        item.flag = 1;
                        gen2drv.put(item);
                        @(gen_ack);
                        break;
                    end
                    else if (i > waiting_time*5 && (i % waiting_time != 0)) begin
                        item.flag = 1;
                        gen2drv.put(item);
                        @(gen_ack);
                        break;
                    end
                    gen2drv.put(item);
                    @(gen_ack);
                    inner1 = 0;
                    inner2 = 1;
                    inner3 = 0;
                end
            end

            //testing flag while edge of wait_timer interval
            for ( i = i+1; i < 2*(waiting_time*10); i++) begin
                item = new();
                if (i > 2*(waiting_time*5) && (i+1) % waiting_time == 0) begin
                    item.flag = 1;
                    gen2drv.put(item);
                    @(gen_ack);
                    $display("i before break:", i);
                    break;
                end
                gen2drv.put(item);
                @(gen_ack);
                inner1 = 0;
                inner2 = 1;
                inner3 = 0;
            end


            item = new();
            item.rst_n      = 1;
            item.start      = 0; 
            item.flag       = 1;
            item.wait_timer = waiting_time;
            inner1 = 0;
            inner2 = 0;
            inner3 = 0;
            gen2drv.put(item);
            @(gen_ack);

            item = new();
            item.rst_n      = 1;
            item.start      = 1; // start signal rising edge
            item.flag       = 1;
            item.wait_timer = waiting_time;
            waiting_time    = waiting_time;
            gen2drv.put(item);
            @(gen_ack);

            // testing max value for count
            for ( i = 0; i < (waiting_time*32 + waiting_time*1); i++) begin
                item = new();
                gen2drv.put(item);
                @(gen_ack);
                inner1 = 0;
                inner2 = 0;
                inner3 = 1;
            end
        end
    endtask

endclass
