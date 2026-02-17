/*`include "env.sv"
`include "and_sequence.sv"
`include "or_sequence.sv"
`include "xor_sequence.sv"
`include "addition_sequence.sv"
`include "seq_item.sv"

class test extends uvm_test;
  `uvm_component_utils(test)
  env env_;
  alu_config cfg;
  virtual alu_if vif;
  addition_sequence add_seq;
  xor_sequence xor_seq;
  and_sequence and_seq;
  or_sequence or_seq;
  function new(string name = "test", uvm_component parent = null);
    super.new(name, parent);
  endfunction
  function void build_phase(uvm_phase phase);
    super.build_phase(phase);
    env_ = env::type_id::create("env_", this);
    cfg =alu_config::type_id::create("cfg",this);
    add_seq = addition_sequence::type_id::create("add_seq", this);
    xor_seq = xor_sequence::type_id::create("xor_seq", this);
    and_seq = and_sequence::type_id::create("and_seq", this);
    or_seq = or_sequence::type_id::create("or_seq", this);
    if (!uvm_config_db #(virtual alu_if)::get(this,"","ALU_IF",cfg.vif))
        `uvm_fatal("NOVIF", "Virtual interface not found in config DB")
        uvm_config_db #(alu_config)::set(this, "*", "vif", cfg);
  endfunction
  task  run_phase (uvm_phase phase);
  super.run_phase(phase);
  phase.raise_objection(this);
  `uvm_info("run_phase","Starting addition test sequences", UVM_LOW)
    add_seq.start(env_.agent_.sequencer_);
     `uvm_info("run_phase","Starting xor test sequences", UVM_LOW)
    xor_seq.start(env_.agent_.sequencer_);
     `uvm_info("run_phase","Starting and test sequences", UVM_LOW)
    and_seq.start(env_.agent_.sequencer_);
     `uvm_info("run_phase","Starting or test sequences", UVM_LOW)
    or_seq.start(env_.agent_.sequencer_);
    phase.drop_objection(this);

    
  endtask //
endclass*/
`include "env.sv"
`include "and_sequence.sv"
`include "or_sequence.sv"
`include "xor_sequence.sv"
`include "addition_sequence.sv"
`include "seq_item.sv"

// ----------------------------------------------------------------------------
// BASE TEST
// Handles configuration, environment creation, and interface passing.
// Does NOT run any sequences.
// ----------------------------------------------------------------------------
class base_test extends uvm_test;
  `uvm_component_utils(base_test)

  env env_;
  alu_config cfg;
  
  function new(string name = "base_test", uvm_component parent = null);
    super.new(name, parent);
  endfunction

  function void build_phase(uvm_phase phase);
    super.build_phase(phase);
    
    // 1. Create Environment and Config
    env_ = env::type_id::create("env_", this);
    cfg  = alu_config::type_id::create("cfg", this);
    
    // 2. Get Interface from Top (matches your previous code)
    if (!uvm_config_db #(virtual alu_if)::get(this, "", "ALU_IF", cfg.vif)) begin
       `uvm_fatal("NOVIF", "Virtual interface not found in config DB")
    end

    // 3. Set Config for Agent
    uvm_config_db #(alu_config)::set(this, "*", "vif", cfg);
  endfunction
  
  // NOTE: run_phase is empty here! 
  // The child classes will fill it in.
endclass


// ----------------------------------------------------------------------------
// TEST 1: ADDITION
// ----------------------------------------------------------------------------
class test_addition extends base_test;
  `uvm_component_utils(test_addition)

  addition_sequence seq;

  function new(string name="test_addition", uvm_component parent=null);
    super.new(name, parent);
  endfunction

  task run_phase(uvm_phase phase);
    seq = addition_sequence::type_id::create("seq");

    phase.raise_objection(this);
    `uvm_info("TEST", "Running ADDITION Sequence", UVM_LOW)
    seq.start(env_.agent_.sequencer_);
    phase.drop_objection(this);
  endtask
endclass


// ----------------------------------------------------------------------------
// TEST 2: XOR
// ----------------------------------------------------------------------------
class test_xor extends base_test;
  `uvm_component_utils(test_xor)

  xor_sequence seq;

  function new(string name="test_xor", uvm_component parent=null);
    super.new(name, parent);
  endfunction

  task run_phase(uvm_phase phase);
    seq = xor_sequence::type_id::create("seq");

    phase.raise_objection(this);
    `uvm_info("TEST", "Running XOR Sequence", UVM_LOW)
    seq.start(env_.agent_.sequencer_);
    phase.drop_objection(this);
  endtask
endclass


// ----------------------------------------------------------------------------
// TEST 3: AND
// ----------------------------------------------------------------------------
class test_and extends base_test;
  `uvm_component_utils(test_and)

  and_sequence seq;

  function new(string name="test_and", uvm_component parent=null);
    super.new(name, parent);
  endfunction

  task run_phase(uvm_phase phase);
    seq = and_sequence::type_id::create("seq");

    phase.raise_objection(this);
    `uvm_info("TEST", "Running AND Sequence", UVM_LOW)
    seq.start(env_.agent_.sequencer_);
    phase.drop_objection(this);
  endtask
endclass


// ----------------------------------------------------------------------------
// TEST 4: OR
// ----------------------------------------------------------------------------
class test_or extends base_test;
  `uvm_component_utils(test_or)

  or_sequence seq;

  function new(string name="test_or", uvm_component parent=null);
    super.new(name, parent);
  endfunction

  task run_phase(uvm_phase phase);
    seq = or_sequence::type_id::create("seq");

    phase.raise_objection(this);
    `uvm_info("TEST", "Running OR Sequence", UVM_LOW)
    seq.start(env_.agent_.sequencer_);
    phase.drop_objection(this);
  endtask
endclass
// ----------------------------------------------------------------------------
// TEST: ALL SEQUENCES (Master Test)
// Runs Addition -> XOR -> AND -> OR in a single simulation
// ----------------------------------------------------------------------------
class test_all extends base_test;
  `uvm_component_utils(test_all)

  // Declare handles for all sequences
  addition_sequence add_seq;
  xor_sequence      xor_seq;
  and_sequence      and_seq;
  or_sequence       or_seq;

  function new(string name="test_all", uvm_component parent=null);
    super.new(name, parent);
  endfunction

  task run_phase(uvm_phase phase);
    // 1. Create all sequence objects
    add_seq = addition_sequence::type_id::create("add_seq");
    xor_seq = xor_sequence::type_id::create("xor_seq");
    and_seq = and_sequence::type_id::create("and_seq");
    or_seq  = or_sequence::type_id::create("or_seq");

    phase.raise_objection(this);

    // 2. Execute them one by one
    `uvm_info("TEST", ">>> Starting ADDITION Sequence <<<", UVM_LOW)
    add_seq.start(env_.agent_.sequencer_);

    `uvm_info("TEST", ">>> Starting XOR Sequence <<<", UVM_LOW)
    xor_seq.start(env_.agent_.sequencer_);

    `uvm_info("TEST", ">>> Starting AND Sequence <<<", UVM_LOW)
    and_seq.start(env_.agent_.sequencer_);

    `uvm_info("TEST", ">>> Starting OR Sequence <<<", UVM_LOW)
    or_seq.start(env_.agent_.sequencer_);

    phase.drop_objection(this);
  endtask
endclass