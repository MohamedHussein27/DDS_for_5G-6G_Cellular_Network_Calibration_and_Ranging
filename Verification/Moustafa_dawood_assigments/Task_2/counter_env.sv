`timescale 1ns/1ps
package counter_env_pkg;

import counter_agent_pkg::*;
import counter_scoreboard_pkg::*;
import counter_coverage_pkg::*;
import counter_seq_item_pkg::*;

class counter_env;

  // Components
  counter_agent      agent;
  counter_scoreboard sb;
  counter_coverage   cov;

  // -----------------------------------
  // Constructor: build only
  // -----------------------------------
  function new();
    agent = new();
    sb    = new();
    cov   = new();
  endfunction

  // -----------------------------------
  // Manual connections (like UVM connect)
  // -----------------------------------
  function void connect(virtual counter_bfm bfm);

    // Agent connections
    agent.connect(bfm);

    // Scoreboard gets data from monitor
    sb.mon2sb_mb = agent.m2sb;

    // Coverage: 2nd approach (monitor-based)
    cov.m2cov_mb = agent.m2cov_mb;

  endfunction

  // -----------------------------------
  // Run all components
  // -----------------------------------
  task run();
    fork
      agent.run();
      sb.run();
      cov.run();
    join_none
  endtask

endclass

endpackage
