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
  

    function new(virtual alu_if vif);
        this.vif = vif;
         // 1. Create Mailboxes
    gen2driv = new();
    mon2scb  = new();
     // 3. Create Components
    gen  = new(gen2driv);
    driv = new(vif, gen2driv);
    cov = new(); 
    mon  = new(vif, mon2scb, cov); 
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