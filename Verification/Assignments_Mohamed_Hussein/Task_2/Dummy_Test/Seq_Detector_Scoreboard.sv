`include "Seq_Detector_State_Classes.sv"
import shared_pkg::*;
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
            // We calculate where to go *next* cycle based on *current* input & current state (mealy FSM)
            next_state = current_state.transition(vif.rst_n, vif.data_in);
            if (~vif.rst_n) begin 
                current_state = next_state; // Immediate transition on reset (asynchronous reset) 
                cs = ns; // update current state variable for checking
            end
            @(posedge vif.clk);
            // state memory
            current_state = next_state;
            cs = ns;


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