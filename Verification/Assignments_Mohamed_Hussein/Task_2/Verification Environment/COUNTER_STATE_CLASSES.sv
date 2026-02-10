import shared_pkg::*;
typedef class State;
typedef class IdleState;
typedef class CountingState;

// classes file variables are defined in shared package (to be displayed)

class State;

    virtual function State transition(logic rst_n, logic start, logic [15:0] wait_timer, logic flag);
        return this; 
    endfunction
    

    virtual function fsm_output_t get_output();
        fsm_output.busy = 0;
        fsm_output.count_value = 5'b00000;
        return fsm_output; // Default output is 0
    endfunction

    virtual function string display_name(); 
        return "Base"; 
    endfunction
endclass

// child classes

class IdleState extends State;
    function string display_name();
        return "IDLE"; 
    endfunction

    function fsm_output_t get_output(); 
        fsm_output.busy = busy_ref;
        fsm_output.count_value = count_ref;
        return fsm_output; 
    endfunction
    
    function State transition(logic rst_n, logic start, logic [15:0] wait_timer, logic flag);
        State next_state; 
        
        if (rst_n == 0) begin
            IdleState idle   = new();
            next_state       = idle;
            ns               = IDLE; // defined in shared package

            busy_ref         = 0;
            count_ref        = 0;
            internal_counter = 0;
        end else if (start) begin
            CountingState counting = new(); 
            next_state = counting;
            ns = COUNTING;

            busy_ref         = 0;
            deassert_count   = 1; // flag for timing matching regarding ref count value (delaying one clock)
            internal_counter = 0;    
        end else begin
            ns               = IDLE;
            next_state       = this;

            busy_ref         = 0;
            //count_ref        = 0;
            internal_counter = 0;    
        end
        
        return next_state;          
    endfunction
endclass


class CountingState extends State;
    function string display_name(); 
        return "COUNTING"; 
    endfunction
    
    // Override output to 1 only for this state
    function fsm_output_t get_output(); 
        fsm_output.busy = busy_ref;
        fsm_output.count_value = count_ref;
        return fsm_output; 
    endfunction

    function State transition(logic rst_n, logic start, logic [15:0] wait_timer, logic flag);
        State next_state;
        
        if (rst_n == 0) begin
            IdleState idle = new();
            next_state = idle;
            ns = IDLE;
            
            busy_ref         = 0;
            count_ref        = 0;
            internal_counter = 0;
        end 
        // if (internal counter is a multiple of wait_timer) & flag is high, or count_value is max, transition to IDLE
        else if ((internal_counter % wait_timer == 0) && (flag == 1 || (fsm_output.count_value == 31 && !deassert_count))) begin
            IdleState idle   = new();
            next_state       = idle;
            ns               = IDLE;
        end 
        else begin
            CountingState counting = new(); 
            next_state             = counting;
            ns                     = COUNTING;
            busy_ref               = 1;
            if (deassert_count) begin
                count_ref          = 0;
                deassert_count     = 0;
            end
            // i just putted this condition before incerementing the counter to avoid incrementing it at the beginning of the state, mimicking NBA behavouir in RTL
            if (internal_counter != 0 && internal_counter % wait_timer == 0) begin
                count_ref++; // Increment count value at each wait_timer interval
            end
            internal_counter++; // Increment internal counter
        end
        
        return next_state;
    endfunction
endclass
