typedef class IdleState;
typedef class wait_state;
typedef enum  {
    IDLE   = 0,
    WAIT     = 1,
    COUNT   = 2
} STATES ;


class State;
    static logic prev_start;
    static logic busy;
    static logic [4:0] count_value=0;
    static logic [15:0] internal_count=0;
    static STATES  current_state;
    static STATES  next_state;
    virtual function State transition(logic start, logic [15:0] wait_timer, logic flag, logic rst_n);
    prev_start = start;
        return this; // Default to stay in the same state
    endfunction

    function State reset_triggered();
        IdleState next_state;
        $display("reset state triggered !!");
        next_state = new();
        return next_state;
    endfunction
    
    function State go_to_idle_state();
        IdleState next_state;
        next_state = new();
        return next_state;
    endfunction

    function State go_to_wait_state();
        wait_state next_state;
        next_state = new();
        return next_state;
    endfunction


endclass : State

class IdleState extends State;
    function State transition(logic start, logic [15:0] wait_timer, logic flag, logic rst_n);
        current_state = IDLE;
        busy = 0;
       // count_value = 0;
        // 1. Check for Reset
        if (!rst_n) begin
            prev_start = 0; // Reset our static history
            next_state = IDLE;
            internal_count = 0;
            count_value = 0;
            return reset_triggered();
        end 
        
        // 2. DETECT RISING EDGE
        else if (start && !prev_start) begin 
            $display("RISING EDGE DETECTED! Going to WAIT.");
            
            
            prev_start = start; 
            
            next_state = WAIT;
            return go_to_wait_state();
        end 
        // 3. Stay in IDLE
        else begin
           
            prev_start = start; 
            
            next_state = IDLE;
            return this;
        end
    endfunction
endclass 
class wait_state extends State;
    function State transition(logic start, logic [15:0] wait_timer, logic flag, logic rst_n);
        current_state = WAIT;
        prev_start = start; // Keep the history up to date even when busy!
        busy = 1;
        $display("the value of start: %0b wait_timer: %0b flag: %0b t:%0p",start,wait_timer,flag,$time);
        if (!rst_n) begin
            count_value = 0;
             internal_count = 0;
            next_state = IDLE;
            return reset_triggered();
        end  else if (wait_timer==internal_count) begin 
            if (!flag && count_value !=31) begin 
            $display("we are in count_Stage and going to wait_state again!");
            next_state = WAIT; 
            internal_count = 1;
            count_value = count_value + 1;
            return go_to_wait_state();
        end 
         else if (flag || count_value ==31) begin 
            $display("we should be in count_Stage and but we are going to idle_state !!");
            next_state = IDLE; 
                internal_count = 0;
            return go_to_idle_state();
        end 
        end 
            else if (wait_timer > internal_count) begin 
                internal_count = internal_count + 1;
        $display("we are in wait_State and still be there");
        next_state = WAIT;
        return this;
        end
        else begin
            $display("we are in wait_State  !!");
            return this;
            end
    endfunction
endclass

