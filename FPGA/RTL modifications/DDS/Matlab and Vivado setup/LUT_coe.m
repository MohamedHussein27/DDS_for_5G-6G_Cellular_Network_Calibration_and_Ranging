% 1. Open the file for reading strings
fid_in = fopen('dds_lut.hex', 'r');
if fid_in == -1
    error('Could not open dds_lut.hex. Ensure the file is in the current folder.');
end

% 2. Use textscan to read the data as strings (%s)
% This prevents MATLAB from trying to interpret "0D" as a decimal number
dataArray = textscan(fid_in, '%s');
fclose(fid_in);

% Extract the cell array of hex strings
hex_data = dataArray{1}; 

% 3. Open the output .coe file
fid_out = fopen('dds_golden_hex.coe', 'w');

% 4. Write Vivado mandatory headers for Hexadecimal
fprintf(fid_out, 'memory_initialization_radix=16;\n');
fprintf(fid_out, 'memory_initialization_vector=\n');

% 5. Write the hex strings separated by commas
for i = 1:(length(hex_data)-1)
    fprintf(fid_out, '%s,\n', hex_data{i});
end

% 6. Write the very last value with a semicolon
fprintf(fid_out, '%s;\n', hex_data{end});

fclose(fid_out);
disp('Successfully created dds_golden_hex.coe using hexadecimal strings!');