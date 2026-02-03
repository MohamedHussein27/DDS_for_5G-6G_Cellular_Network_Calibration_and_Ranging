`include "alu_item.sv"
class alu_generator;
  mailbox #(alu_item) gen2driv;
  function new(mailbox gen2driv);
    this.gen2driv = gen2driv;
  endfunction

  task run();
    alu_item item;
    for (int op = 0; op < 4; op++) begin
      for (int a = 0; a < 16; a++) begin
        for (int b = 0; b < 16; b++) begin
            
          item = new();
          item.op = operation_t'(op); 
          item.a = a;
          item.b = b;
          gen2driv.put(item);
          
        end
      end
    end
  endtask
  
endclass