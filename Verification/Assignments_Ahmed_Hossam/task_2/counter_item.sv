`ifndef COUNTER_ITEM_SV
`define COUNTER_ITEM_SV   // <--- Fixed typo (M instead of N)

class counter_item;       // <--- Renamed to match your Scoreboard
    
    // ------------------------------------------------
    // 1. INPUTS (Randomized for Driver)
    // ------------------------------------------------
    rand bit        rst_n;
    rand bit        start;
    //rand bit [15:0] wait_timer;
    bit [15:0] wait_timer=16;
    rand bit        flag;

    // ------------------------------------------------
    // 2. OUTPUTS (Captured by Monitor for Scoreboard)
    // ------------------------------------------------
    // These are NOT 'rand' because they are results from the DUT.
    bit             busy;
    bit [4:0]       count_value;

    // ------------------------------------------------
    // 3. CONSTRAINTS
    // ------------------------------------------------
    constraint cons {
        start dist {0:=40, 1:=60};
        rst_n dist {0:=1,  1:=99};
        
        // Optional: Constrain timer to smaller values during debug 
        // so you don't wait 65,000 cycles for one test!
        //wait_timer inside {[1:20]}; 
    }
    
    // Optional: A helper function to print meaningful logs
    function string convert2string();
        return $sformatf("rst_n=%b start=%b wait=%0d flag=%b | busy=%b count=%0d", 
                         rst_n, start, wait_timer, flag, busy, count_value);
    endfunction

endclass
`endif