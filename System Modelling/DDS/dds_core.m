function dds_out = dds_core(M, Nacc, LUT_bits, DT_Mode) %#codegen
%DDS_CORE Hardware-accurate DDS model

    % -----------------------------------------------------
    % 1. Load Types
    % -----------------------------------------------------
    T = mytypes(DT_Mode);

    assert(LUT_bits <= Nacc, 'LUT_bits must be <= Nacc');

    % -----------------------------------------------------
    % 2. Cast Inputs
    % -----------------------------------------------------
    M_cast = cast(M, 'like', T.M);
    Ns = length(M);

    % -----------------------------------------------------
    % 3. Sine LUT (Phase-to-Amplitude Converter)
    % -----------------------------------------------------
    LUT_size = 2^LUT_bits;
    addr_ROM = 0:LUT_size-1;
    sine_LUT = cast(sin(2*pi*addr_ROM/LUT_size), 'like', T.sine_LUT);

    % -----------------------------------------------------
    % 4. Phase Accumulator
    % -----------------------------------------------------
    acc   = uint64(0);
    phase = zeros(1, Ns, 'uint64');

    for n = 1:Ns
        % Read current phase
        phase(n) = acc;

        % Update accumulator (wrap naturally)
        acc = acc + uint64(M_cast(n));
        acc = bitand(acc, uint64(2^Nacc - 1));
    end

    % -----------------------------------------------------
    % 5. Phase Truncation
    % -----------------------------------------------------
    addr_shift = Nacc - LUT_bits;
    lut_addr = bitshift(phase, -addr_shift);
    lut_addr = cast(lut_addr, 'like', T.lut_addr);

    % -----------------------------------------------------
    % 6. ROM Lookup
    % -----------------------------------------------------
    dds_out = zeros(1, Ns, 'like', T.dds_out);
    dds_out(:) = sine_LUT(double(lut_addr) + 1);
end
