# Open the COE file and read the lines
with open('rtl_rx_in_re.coe', 'r') as f:
    lines = f.readlines()

# Open a new Verilog file to write the output
with open('rx_input_rom_hardcoded.v', 'w') as out:
    # Write the Verilog header
    out.write("`timescale 1ns / 1ps\n\n")
    out.write(
        "module rx_input_rom_hardcoded #(\n    parameter WL = 16,\n    parameter DEPTH = 4096\n)(\n")
    out.write("    input  wire                      clk,\n")
    out.write(
        "    input  wire [$clog2(DEPTH)-1:0]  addr,\n    output reg  signed [WL-1:0]      data_out_re\n);\n\n")
    out.write(
        "    (* ram_style = \"block\" *) reg [WL-1:0] rom_re [0:DEPTH-1];\n\n")
    out.write("    initial begin\n")

    # Extract the hex values and write the 4096 lines
    index = 0
    for line in lines:
        line = line.strip()
        # Ignore the COE header lines
        if line.startswith('memory') or line == '':
            continue

        # Remove the comma or semicolon at the end of the hex value
        hex_val = line.replace(',', '').replace(';', '')

        out.write(f"        rom_re[{index}] = 16'h{hex_val};\n")
        index += 1

    # Write the Verilog footer
    out.write("    end\n\n")
    out.write("    always @(posedge clk) begin\n")
    out.write("        data_out_re <= rom_re[addr];\n")
    out.write("    end\n\nendmodule\n")

print(f"Successfully generated rx_input_rom_hardcoded.v with {index} samples!")

# # Open the COE file and read the lines
# with open('rtl_rx_in_im.coe', 'r') as f:
#     lines = f.readlines()

# # Open a new Verilog file to write the output
# with open('rx_input_rom_hardcoded.v', 'w') as out:
#     # Write the Verilog header
#     out.write("`timescale 1ns / 1ps\n\n")
#     out.write(
#         "module rx_input_rom_hardcoded #(\n    parameter WL = 16,\n    parameter DEPTH = 4096\n)(\n")
#     out.write("    input  wire                      clk,\n")
#     out.write(
#         "    input  wire [$clog2(DEPTH)-1:0]  addr,\n    output reg  signed [WL-1:0]      data_out_im\n);\n\n")
#     out.write(
#         "    (* ram_style = \"block\" *) reg [WL-1:0] rom_im [0:DEPTH-1];\n\n")
#     out.write("    initial begin\n")

#     # Extract the hex values and write the 4096 lines
#     index = 0
#     for line in lines:
#         line = line.strip()
#         # Ignore the COE header lines
#         if line.startswith('memory') or line == '':
#             continue

#         # Remove the comma or semicolon at the end of the hex value
#         hex_val = line.replace(',', '').replace(';', '')

#         out.write(f"        rom_im[{index}] = 16'h{hex_val};\n")
#         index += 1

#     # Write the Verilog footer
#     out.write("    end\n\n")
#     out.write("    always @(posedge clk) begin\n")
#     out.write("        data_out_im <= rom_im[addr];\n")
#     out.write("    end\n\nendmodule\n")

# print(f"Successfully generated rx_input_rom_hardcoded.v with {index} samples!")
