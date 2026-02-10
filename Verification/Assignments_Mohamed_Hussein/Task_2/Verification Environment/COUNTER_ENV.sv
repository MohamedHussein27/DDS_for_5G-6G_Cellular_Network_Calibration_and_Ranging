`include "COUNTER_AGENT.sv"
`include "COUNTER_COVERAGE.sv"
`include "COUNTER_SCOREBOARD.sv"
import shared_pkg::*;
class counter_env;
    string name;

    // virtual interface
    virtual counter_if counter_vif;

    // agent, scoreboard and coverage class objects
    counter_agent agt;
    counter_scoreboard sb;
    counter_coverage cov;

    function new(string name = "counter_env");
        this.name = name;
        agt = new();
        sb  = new();
        cov = new();
    endfunction

    // function to connect seq_item to scoreboard and coverage
    function void connect();
        // connect virtual interface to agent virtual interface
        agt.counter_vif = counter_vif;
        // connect mailbox to scoreboard and coverage
        sb.mon2sb = agt.mon2sb;
        cov.mon2cov = agt.mon2cov;
    endfunction

    task run();
        agt.connect(); // to connect mailboxes and vifs
        fork
            agt.run();
            sb.run();
            cov.run();
        join_any
        sb.report(); // call report task to print the final count of errors and correct transactions
    endtask

endclass