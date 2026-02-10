`include "counter_agent.sv"
`include "counter_scoreboard.sv"
class counter_env ;
counter_scoreboard sb;
counter_agent agent;
virtual counter_if vif;
mailbox #(counter_item) mon2scb;
function new();    
endfunction
function void connecting(virtual counter_if vif);
this.vif=vif;
agent=new();
agent.connecting(vif);
sb=new();//////////////////////////
sb.connecting(agent.mon2scb);
endfunction
function run();
fork
agent.run();
sb.run();
    
join_none

endfunction
endclass