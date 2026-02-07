`include "alu_generator.sv" 
`include "alu_driver.sv"
`include "alu_monitor.sv"
class agent;
    alu_generator  gen;
    alu_driver     driv;
    alu_monitor    mon;
    alu_coverage   cov;  // The coverage component
  
    // Mailboxes & Interface
  mailbox #(alu_item) gen2driv;
  mailbox #(alu_item) mon2scb;
  virtual alu_if vif;
  function new();   
  endfunction

    function void connecting(virtual alu_if vif);
        this.vif = vif;
         // 1. Create Mailboxes
    gen2driv = new();
    mon2scb  = new();
     // 3. Create Components
    gen  = new();
    gen.connecting(gen2driv); // Connect Generator to Driver's mailbox
    driv = new();
    driv.connecting(vif, gen2driv); // Connect Driver to Generator's mailbox
    cov = new(); 
    mon  = new();
    mon.connecting(vif, mon2scb, cov); // Connect Monitor to Interface & Scoreboard's mailbox 
    endfunction //new()

    task run();
    fork
      driv.run();
      mon.run();
      
    join_none 
    
    // 5. Run Generator (Main Thread)
    gen.run(); 
    endtask //run()

endclass //agent