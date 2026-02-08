import shared_pkg::*;
typedef class State;
typedef class IdleState;
typedef class S1State;
typedef class S10State;
typedef class S101State;

class State;
    // Transition Logic (Next State)
    virtual function State transition(logic data_in);
        return this; 
    endfunction
    
    // Output Logic (New Request)
    virtual function logic get_output();
        return 0; // Default output is 0
    endfunction

    virtual function string get_name(); 
        return "Base"; 
    endfunction
endclass

// child classes

class IdleState extends State;
    function string get_name();
        return "IDLE"; 
    endfunction
    
    function State transition(logic data_in);
        State next_state; 
        
        if (data_in) begin
            S1State s1 = new(); 
            next_state = s1;
            ns = S1; // defined in shared package     
        end else begin
            ns = IDLE;
            next_state = this;      
        end
        
        return next_state;          
    endfunction
endclass

class S1State extends State;
    function string get_name(); 
        return "S1"; 
    endfunction
    
    function State transition(logic data_in);
        State next_state;
        
        if (data_in) begin
            ns = S1;
            next_state = this;      // Stay in S1
        end else begin
            S10State s10 = new();
            next_state = s10;
            ns = S10;
        end
        
        return next_state;
    endfunction
endclass

class S10State extends State;
    function string get_name(); 
        return "S10"; 
    endfunction
    
    function State transition(logic data_in);
        State next_state;
        
        if (data_in) begin
            S101State s101 = new();
            next_state = s101;
            ns = S101;
        end else begin
            IdleState idle = new();
            next_state = idle;
            ns = IDLE;
        end
        
        return next_state;
    endfunction
endclass

class S101State extends State;
    function string get_name(); 
        return "S101"; 
    endfunction
    
    // Override output to 1 only for this state
    function logic get_output(); 
        return 1; 
    endfunction

    function State transition(logic data_in);
        State next_state;
        
        if (data_in) begin
            S1State s1 = new();
            next_state = s1;
            ns = S1;
        end else begin
            IdleState idle = new();
            next_state = idle;
            ns = IDLE;
        end
        
        return next_state;
    endfunction
endclass