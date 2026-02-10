`include "alu_scoreboard.sv"
`include "agent.sv"
class alu_env;
    // Components
  agent         agent_; 
  alu_scoreboard scb;
  

  mailbox #(alu_item) mon2scb;
  virtual alu_if vif;
  
  // Constructor
  function new();   
  endfunction
  function void connecting(virtual alu_if vif);
    this.vif = vif;
    mon2scb  = new();
    
    // 2. Create Coverage FIRST (So we can pass it to Monitor)
    
    agent_ = new();
    agent_.connecting(vif); // Connect Agent's components   
    scb  = new();
    scb.connecting(agent_.mon2scb); // Connect Scoreboard to Monitor's mailbox
  endfunction
  
  // Run Task
  task run();
    fork
      agent_.run();
      scb.run();
    join_none // changed to join_none so we proceed to gen.run()
    
  
    
    
  endtask
  
endclass