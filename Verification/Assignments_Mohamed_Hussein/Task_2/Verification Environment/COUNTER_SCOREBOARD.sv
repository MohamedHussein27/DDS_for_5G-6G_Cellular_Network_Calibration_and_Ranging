import shared_pkg::*;
import counter_seq_item_pkg::*;
`include "COUNTER_STATE_CLASSES.sv"
class counter_scoreboard;
    string name;

    mailbox #(counter_seq_item) mon2sb;

    // State Handles
    State current_state; 
    State next_state;

    function new(string name = "SCOREBOARD");
        this.name = name;
        mon2sb = new();
    endfunction

    // compare function
    function void check_data (counter_seq_item tr, fsm_output_t fsm_out);
        if (tr.count_value !== fsm_out.count_value) begin
            error_count_out++;
            $display("error in fsm output, ref_count_value is: %0d     while dut count value is: %0d", fsm_out.count_value, tr.count_value  );
        end
        else
            correct_count_out++;
        if(tr.busy !== fsm_out.busy) begin
            error_count_busy++;
            $display("error in fsm output, ref_busy is: %0d    while dut busy is: %0d", fsm_out.busy, tr.busy  );
        end
        else
            correct_count_busy++;
    endfunction


    task run();
        // Initialize
        counter_seq_item item;
        
        IdleState idle = new();
        current_state = idle;
        next_state = current_state;

        forever begin
            item = new();
            mon2sb.get(item);
            //$display("we don't get anything" , item.start);

            // We calculate where to go *next* cycle based on *current* input & current state (mealy FSM)
            next_state = current_state.transition(item.rst_n, item.start, item.wait_timer, item.flag);
            
            // state memory
            current_state = next_state;
            cs = ns;

            // calling get_output fn
            fsm_out = current_state.get_output(); // fsm_out is defined in shared package

            check_data(item, fsm_out);
        end
    endtask

    task report(); // this task is the counter of report phase in UVM
        $display("error count out = %0d, correct count out = %0d", error_count_out, correct_count_out);
        $display("error count busy = %0d, correct count busy = %0d", error_count_busy, correct_count_busy);
    endtask
endclass
