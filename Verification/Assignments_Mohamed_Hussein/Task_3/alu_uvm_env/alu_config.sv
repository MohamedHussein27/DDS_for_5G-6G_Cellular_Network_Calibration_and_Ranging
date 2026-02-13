package alu_config_obj_pkg;
    import uvm_pkg::*;
    `include "uvm_macros.svh"
    class alu_config_obj extends uvm_object;
        `uvm_object_utils(alu_config_obj)

        // virtual interface (main)
        virtual alu_if alu_vif;

        // constructor
        function new(string name = "alu_config_obj");
            super.new(name);
        endfunction
    endclass
endpackage
