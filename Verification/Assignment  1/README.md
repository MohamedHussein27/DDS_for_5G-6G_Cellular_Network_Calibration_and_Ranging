# ALU Verification Ramp-Up Project

## Overview
This repository contains a complete functional verification environment for a combinational Arithmetic Logic Unit (ALU).  
The ALU design is intentionally injected with a large number of functional bugs to simulate an industrial-style verification ramp-up scenario.

The goal of this project is to demonstrate systematic bug discovery, verification methodology, and coverage-driven validation using SystemVerilog.

---

## Design Under Test (DUT)
The DUT is a purely combinational ALU that supports basic arithmetic and logical operations.  
Due to the injected bugs, the DUT exhibits multiple functional and structural mismatches relative to the expected ALU specification, making it suitable for verification-focused evaluation rather than design correctness.

---

## Verification Methodology
•	Verification approach: (random + constrained-random)
•	In my verification flow I found some mismatches between golden model and dut, but for every mismatch with the same op the out signal is something unique, so I decided to take five mismatch cases with the same opcode for each operation (ADD, XOR, AND, OR), comparing the inputs and outputs trying to find a relation between them to catch the bugs.

No clocking mechanism is used, as the DUT is purely combinational. Instead, temporal separation between stimulus driving and sampling is enforced to avoid race conditions.

---

## Verification Environment Architecture

> **Verification Environment Block Diagram**  

![Verification Environment Architecture](https://github.com/MohamedHussein27/Verification_Environment_for_an_ALU/blob/main/Documentation/Env_Structure.png)

---


## Key Verification Components
- **Driver:** Generates and applies stimulus transactions to the DUT
- **Monitor:** Observes DUT inputs and outputs and converts them into transactions
- **Reference Model:** Computes expected ALU results based on the specification
- **Scoreboard:** Compares DUT outputs against reference model predictions
- **Coverage Collector:** Measures functional coverage to assess verification completeness
- **Sequence Item:** Defines a single ALU transaction shared across verification components

---

## Results
The verification environment successfully detects a large number of functional bugs across different ALU operations.  
Coverage results and scoreboard statistics are used to evaluate verification completeness and effectiveness.


---

## Notes
This project focuses on verification quality and methodology rather than DUT correctness.  
The DUT behavior should not be considered a reference for a correct ALU implementation.

---

## Author
Mohamed Hussein  
Graduation Project – Verification Track

