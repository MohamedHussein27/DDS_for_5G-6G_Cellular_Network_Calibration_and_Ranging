package alu_config_pkg;
	import uvm_pkg::*;
    import alu_pkg::*;
	`include "uvm_macros.svh"

	class alu_config extends uvm_object;
		`uvm_object_utils(alu_config)

		virtual alu_if alu_vif;

		function new(string name = "alu_config");
			super.new(name);
		endfunction 
		
	endclass 
	
endpackage 