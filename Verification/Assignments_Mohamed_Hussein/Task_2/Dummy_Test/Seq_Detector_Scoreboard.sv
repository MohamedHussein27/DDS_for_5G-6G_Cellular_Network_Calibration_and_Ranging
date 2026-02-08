`include "Seq_Detector_State_Classes.sv"
class my_scoreboard;
    
    // Virtual Interface to access signals
    virtual seq_intf vif;

    // State Handles
    State current_state; 
    State next_state;

    function void set_vif(virtual seq_intf vif);
        this.vif = vif;
    endfunction


    task run();
        // Initialize
        IdleState idle = new();
        current_state = idle;
        next_state = current_state;

        forever begin
            @(posedge vif.clk);
            

            if (!vif.rst_n) begin
                // Reset Logic
                IdleState idle = new();
                next_state = idle;
                current_state = next_state;
            end 
            else begin
                // state memory
                // We calculate where to go *next* cycle based on *current* input & current state (mealy FSM)
                next_state = current_state.transition(vif.data_in); 
                current_state = next_state;
            end

            // checking
            @(negedge vif.clk); // Check output at the end of the cycle
            if (vif.seq_detected !== current_state.get_output()) begin
                $error("Mismatch! State: %s | DUT: %b | Model: %b", 
                        current_state.get_name(), vif.seq_detected, current_state.get_output());
            end else begin
                $display("Time: %0t | State: %s | Match: OK", $time, current_state.get_name());
            end
        end
    endtask
endclass