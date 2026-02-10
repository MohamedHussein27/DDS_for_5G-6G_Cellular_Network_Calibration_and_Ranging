`include "counter_item.sv"
// Include your FSM class file here so the scoreboard recognizes 'State', 'IdleState', etc.
`include "golden_model.sv" 

class counter_scoreboard;
    mailbox #(counter_item) mon2scb;
    int error_count = 0;
    int pass_count = 0;
    
    // 1. Declare the Golden Model Handle
    // This handle must persist across every 'check' call to remember the state.
    State fsm_model; 

    function new();
        // 2. Initialize the Model (Start in IDLE)
        IdleState init_state = new();
        fsm_model = init_state;
    endfunction

    function void connecting (mailbox mon2scb);
        this.mon2scb = mon2scb;
    endfunction

    task run();
        forever begin
            counter_item item;
            mon2scb.get(item);
//////
if (item.rst_n==0) begin
    fsm_model = fsm_model.reset_triggered();
end


//////
            check(item);
        end
    endtask

    task check(counter_item item);
        bit [4:0] expected_count_value;
        bit expected_busy_flag;

        // ============================================================
        // 3. GOLDEN MODEL LOGIC
        // ============================================================
        // We feed the monitor's data into the model to step it forward one clock cycle.
        // The model returns the new state object (Idle, Wait, or Count).
        
        fsm_model = fsm_model.transition(
            item.start, 
            item.wait_timer, 
            item.flag, 
            item.rst_n
        );
       // #1;
        // 4. Extract Expected Outputs from the Model
        expected_count_value = fsm_model.count_value;
        expected_busy_flag   = fsm_model.busy;

        // ============================================================
        // 5. COMPARISON LOGIC
        // ============================================================
        
        if (item.count_value != expected_count_value || item.busy != expected_busy_flag) begin
            error_count++;
            $display("[SCOREBOARD] ERROR @%0t", $time);
            $display("    Inputs  -> Start: %b | Wait: %0d | Flag: %b", item.start, item.wait_timer, item.flag);
            $display("    Expected-> Busy: %b | Count: %0d", expected_busy_flag, expected_count_value);
            $display("    Actual  -> Busy: %b | Count: %0d", item.busy, item.count_value);
            $display("    State   -> %p", fsm_model.current_state); // Optional: Print current state enum
            $display("    Total Errors So Far: %0d", error_count);
        end else begin
            pass_count++;
            // Optional: Only print every 100 passes to avoid cluttering the log
             $display("[SCOREBOARD] PASS @%0t | Count: %0d", $time, pass_count);
        end
    endtask
endclass