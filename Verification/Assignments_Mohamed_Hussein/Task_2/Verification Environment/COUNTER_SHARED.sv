// this package has the shared signals among packages
package shared_pkg;
    // counters
    int error_count_out     = 0; 
    int correct_count_out   = 0;
    int error_count_busy    = 0;
    int correct_count_busy  = 0;
    // enum
    typedef enum  {
            IDLE          ,
            START         ,
            COUNTING      
        } state_e; 

    typedef struct packed {
        logic busy;
        logic [4:0] count_value; 
    } fsm_output_t;

    // classes file variables
    state_e cs, ns;
    fsm_output_t fsm_output, fsm_out;

    // global variable in classes file
    int internal_counter = 0;
    int count_ref        = 0;
    bit busy_ref         = 0;
    bit deassert_count   = 0; 
endpackage