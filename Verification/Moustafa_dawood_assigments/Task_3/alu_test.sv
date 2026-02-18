package alu_test_pkg;

  import uvm_pkg::*;
  import alu_config_pkg::*;
  import alu_env_pkg::*;
  import or_sequence_pkg::*;
  import and_sequence_pkg::*;
  import xor_sequence_pkg::*;
  import add_sequence_pkg::*;
  import random_op_pkg::*;
  import rst_sequence_pkg::*;

  `include "uvm_macros.svh"

  // ==========================================================
  // BASE TEST
  // ==========================================================
  class base_test extends uvm_test;

    `uvm_component_utils(base_test)

    alu_env    env;
    alu_config alu_cfg;

    function new(string name="base_test", uvm_component parent=null);
      super.new(name,parent);
    endfunction

    function void build_phase(uvm_phase phase);
      super.build_phase(phase);

      env     = alu_env::type_id::create("env", this);
      alu_cfg = alu_config::type_id::create("alu_cfg");

      if(!uvm_config_db#(virtual alu_if)::get(this, "", "ALU_IF", alu_cfg.alu_vif))
        `uvm_fatal("BUILD", "Unable to get ALU interface")

      uvm_config_db#(alu_config)::set(this, "*", "CFG", alu_cfg);

    endfunction

  endclass


  // ==========================================================
  // RESET TEST
  // ==========================================================
  class rst_test extends base_test;

    `uvm_component_utils(rst_test)

    // ✅ ADDED CONSTRUCTOR
    function new(string name="rst_test", uvm_component parent=null);
      super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
      rst_sequence seq;
      phase.raise_objection(this);
      seq = rst_sequence::type_id::create("seq");
      seq.start(env.agt.sqr);
      phase.drop_objection(this);
    endtask

  endclass


  // ==========================================================
  // XOR TEST
  // ==========================================================
  class xor_test extends base_test;

    `uvm_component_utils(xor_test)

    // ✅ ADDED CONSTRUCTOR
    function new(string name="xor_test", uvm_component parent=null);
      super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
      xor_sequence seq;
      phase.raise_objection(this);
      seq = xor_sequence::type_id::create("seq");
      seq.start(env.agt.sqr);
      phase.drop_objection(this);
    endtask

  endclass


  // ==========================================================
  // AND TEST
  // ==========================================================
  class and_test extends base_test;

    `uvm_component_utils(and_test)

    // ✅ ADDED CONSTRUCTOR
    function new(string name="and_test", uvm_component parent=null);
      super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
      and_sequence seq;
      phase.raise_objection(this);
      seq = and_sequence::type_id::create("seq");
      seq.start(env.agt.sqr);
      phase.drop_objection(this);
    endtask

  endclass


  // ==========================================================
  // OR TEST
  // ==========================================================
  class or_test extends base_test;

    `uvm_component_utils(or_test)

    // ✅ ADDED CONSTRUCTOR
    function new(string name="or_test", uvm_component parent=null);
      super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
      or_sequence seq;
      phase.raise_objection(this);
      seq = or_sequence::type_id::create("seq");
      seq.start(env.agt.sqr);
      phase.drop_objection(this);
    endtask

  endclass


  // ==========================================================
  // ADD TEST
  // ==========================================================
  class add_test extends base_test;

    `uvm_component_utils(add_test)

    // ✅ ADDED CONSTRUCTOR
    function new(string name="add_test", uvm_component parent=null);
      super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
      add_sequence seq;
      phase.raise_objection(this);
      seq = add_sequence::type_id::create("seq");
      seq.start(env.agt.sqr);
      phase.drop_objection(this);
    endtask

  endclass


  // ==========================================================
  // RANDOM TEST
  // ==========================================================
  class random_test extends base_test;

    `uvm_component_utils(random_test)

    // ✅ ADDED CONSTRUCTOR
    function new(string name="random_test", uvm_component parent=null);
      super.new(name, parent);
    endfunction

    task run_phase(uvm_phase phase);
      random_op seq;
      phase.raise_objection(this);
      seq = random_op::type_id::create("seq");
      seq.start(env.agt.sqr);
      phase.drop_objection(this);
    endtask

  endclass

endpackage

/*
package alu_test_pkg;
import uvm_pkg::*;
import alu_config_pkg::*;
import alu_env_pkg::*;
import or_sequence_pkg::*;
import and_sequence_pkg::*;
import xor_sequence_pkg::*;
import add_sequence_pkg::*;
import random_op_pkg::*;
import rst_sequence_pkg::*;

`include "uvm_macros.svh"

class alu_test extends uvm_test;
`uvm_component_utils(alu_test)

alu_env env;
alu_config alu_cfg;
xor_sequence xor_seq;
and_sequence and_seq;
or_sequence or_seq;
add_sequence add_seq;
random_op rand_seq;
rst_sequence reset_seq;


    function new(string name = "alu_test ", uvm_component parent = null);
        super.new(name, parent);
        
    endfunction //new()

    function  void build_phase(uvm_phase phase);
        super.build_phase(phase);
        `uvm_info("BUILD_PHASE", "Building the test", UVM_LOW)
        env = alu_env::type_id::create("env", this);
        alu_cfg=alu_config::type_id::create("alu_cfg");
        xor_seq = xor_sequence::type_id::create("xor_seq");
        and_seq = and_sequence::type_id::create("and_seq");
        or_seq = or_sequence::type_id::create("or_seq");
        add_seq = add_sequence::type_id::create("add_seq");
        rand_seq = random_op::type_id::create("rand_seq");
        reset_seq = rst_sequence::type_id::create("reset_seq");




         if(!uvm_config_db#(virtual alu_if)::get(this, "", "ALU_IF", alu_cfg.alu_vif)) //get the vif and assign it to the vif of the cfg
				`uvm_fatal("build_phase", "Test - Unable to get the virtual interface of the ALU from the uvm_config_db")

			uvm_config_db#(alu_config)::set(this, "*", "CFG", alu_cfg);

        
        
    endfunction

    task run_phase(uvm_phase phase);
            super.run_phase(phase);
			phase.raise_objection(this);

            // rest sequnce 
            `uvm_info("run_phase", "Reset Asserted", UVM_LOW)
            reset_seq.start(env.agt.sqr);
            `uvm_info("run_phase", "Reset Deasserted", UVM_LOW)


            // xor sequence
            `uvm_info("run_phase", "Starting XOR Sequence", UVM_LOW)
            xor_seq.start(env.agt.sqr);
             `uvm_info("run_phase", "XOR Sequence completed", UVM_LOW)

            // and sequence
            `uvm_info("run_phase", "Starting AND Sequence", UVM_LOW)
            and_seq.start(env.agt.sqr);
            `uvm_info("run_phase", "AND Sequence completed", UVM_LOW)

            // or sequence
            `uvm_info("run_phase", "Starting OR Sequence", UVM_LOW)
            or_seq.start(env.agt.sqr);
             `uvm_info("run_phase", "OR Sequence completed", UVM_LOW)


            // add sequence
            `uvm_info("run_phase", "Starting ADD Sequence", UVM_LOW)
            add_seq.start(env.agt.sqr);
             `uvm_info("run_phase", "ADD Sequence completed", UVM_LOW)

             // random op sequence
                `uvm_info("run_phase", "Starting Random Sequence", UVM_LOW)
            rand_seq.start(env.agt.sqr);
             `uvm_info("run_phase", "Random Sequence completed", UVM_LOW)



            phase.drop_objection(this);
    endtask : run_phase

   



endclass  
    
endpackage
*/
