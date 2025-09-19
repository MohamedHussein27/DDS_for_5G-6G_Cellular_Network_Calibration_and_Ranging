# DDS for 5G/6G ISAC Calibration and Ranging

## Project Overview

This graduate project, **sponsored by Analog Devices Inc. (ADI)**, focuses on the hardware implementation of a **Direct Digital Synthesis (DDS)** waveform generator and correlator. This IP core is a critical component for **Integrated Sensing and Communication (ISAC)** in next-generation 5G and 6G cellular networks.

ISAC enables a single wireless infrastructure to simultaneously perform high-data-rate communication (e.g., voice, video) and high-resolution sensing (e.g., localization, ranging, and environment mapping). This dual functionality is key to enabling advanced use cases like autonomous vehicles, smart cities, and augmented reality.

The project involves the design of a flexible DDS engine capable of generating wideband waveforms (LFM/NLFM for sensing and OFDM for communication) required for calibrating massive MIMO antenna arrays and performing precise ranging operations in modern networks that support bandwidths up to 1600 MHz.

## Project Goals

*   **System Analysis:** To understand the applications and necessary features for joint integrated sensing and communication (ISAC) in 5G/6G networks.
*   **Algorithm Comparison:** To compare different ISAC signaling schemes (OFDM-based, LFM-based, Hybrid) in terms of algorithmic performance and hardware implementation complexity.
*   **FPGA Implementation:** To design and implement selected ISAC functionalities of the low-physical layer and Digital Front-End (DFE), specifically a programmable DDS and correlator, on an FPGA platform.

## Technical Specifications

*   **Target Application:** 5G/6G Cellular Network Calibration and Ranging
*   **Supported Bandwidths:** 100 MHz, 200 MHz, 400 MHz, up to 1600 MHz
*   **Key Features:**
    *   Programmable waveform generation (LFM, NLFM, OFDM)
    *   High-resolution frequency tuning
    *   Integrated digital correlator for processing reflected sensing signals
    *   SPI interface for control and configuration
*   **Design Flow:**
    *   **System-Level Modeling:** MATLAB/C++
    *   **RTL Design:** Verilog/SystemVerilog
    *   **Verification:** **Cocotb** and **PyUVM** for rigorous functional verification
    *   **Implementation Target:** FPGA and/or ASIC Prototyping

## Project Team

* **Team Leader:** Mohamed Hussein
* **Team Members:**
  * Ahmed Hossam
  * Ahmed Haitham
  * Mohamed Adel
  * Abdelrahman Khaled
  * Moustafa Dawood
  * Youssef Mohamed

## Acknowledgments

This project is developed in collaboration with and sponsored by **Analog Devices Inc. (ADI)**.

## Project Status

**Status:** Active Development 🚧
*   [ ] System-level modeling and algorithm selection (MATLAB/C++)
*   [ ] DDS core RTL design (Verilog/SV)
*   [ ] Correlator IP design (Verilog/SV)
*   [ ] SPI Control Interface implementation
*   [ ] Cocotb/PyUVM testbench development
*   [ ] FPGA Implementation & Validation


## References

[1] 3GPP TS 38.104, "Base Station (BS) radio transmission and reception," 2022.
[2] L. L. et al., "A Survey of the Functional Splits Proposed for 5G Mobile Crosshaul Networks," IEEE Communications Surveys & Tutorials, vol. 21, no. 1, 2019.
[3] M. P. et al., "Understanding O-RAN: Architecture, Interfaces, Algorithms, Security, and Research Challenges," ArXiv 2022.
[4] Analog Devices, Inc., "Fundamentals of Direct Digital Synthesis (DDS)," MT-085 Tutorial.
[5] Analog Devices, Inc., "All About Direct Digital Synthesis," Analog Dialogue.
