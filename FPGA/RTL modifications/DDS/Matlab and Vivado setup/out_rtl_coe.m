% 1. Open the decimal file for reading
fid_in = fopen('rtl_dds_out.txt', 'r');
if fid_in == -1
    error('Could not open file. Ensure rtl_dds_out.txt is in the current folder.');
end

% 2. Read data as strings to handle scientific notation safely
dataArray = textscan(fid_in, '%s');
fclose(fid_in);

% 3. Convert strings to numeric values and round to nearest integer
data = str2double(dataArray{1});
data = round(data);

% 4. Apply 8-bit Two's Complement for negative numbers
% This converts -17 to 239, etc.
data(data < 0) = data(data < 0) + 256;

% 5. Create and open the new hex .coe file
fid_out = fopen('rtl_dds_golden_hex.coe', 'w');

% 6. Write Vivado mandatory headers for Hexadecimal (Radix 16)
fprintf(fid_out, 'memory_initialization_radix=16;\n');
fprintf(fid_out, 'memory_initialization_vector=\n');

% 7. Write the data in Hex format (%02X forces 2-digit uppercase hex)
for i = 1:(length(data)-1)
    fprintf(fid_out, '%02X,\n', data(i));
end

% 8. Write the very last value with a semicolon
fprintf(fid_out, '%02X;\n', data(end));

fclose(fid_out);
disp('Successfully created rtl_dds_golden_hex.coe!');