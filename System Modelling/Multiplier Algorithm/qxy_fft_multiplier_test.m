%% ---------------------------------------------------------
% Testbench for Qx.y FFT Multiplier (double, single, fixed)
% ---------------------------------------------------------
clear; clc;

%% Parameters
Npts = 2048;
x = 1;           % Integer bits including sign
y = 10;          % Fractional bits
rng(1);          % Repeatability

%% Generate random FFT-like inputs (bounded to Qx.y)
FS = 0.95 * 2^(x-1);      % Safe full-scale bound

A = (rand(Npts,1)-0.5)*2*FS + 1j*(rand(Npts,1)-0.5)*2*FS;
B = (rand(Npts,1)-0.5)*2*FS + 1j*(rand(Npts,1)-0.5)*2*FS;


%% ----------------------------
% Call hardware-style multiplier
%% ----------------------------
Y_double = qxy_fft_multiplier(A, B, x, y, 'double');
Y_single = qxy_fft_multiplier(A, B, x, y, 'single');
Y_fixed  = qxy_fft_multiplier(A, B, x, y, 'fixed');

%% Reference (double precision)
Y_ref = A .* B;

%% Convert fixed-point to double for error analysis
Y_fixed_dbl = double(Y_fixed);

%% ----------------------------
% Error analysis
%% ----------------------------
err_double = Y_double - Y_ref;
err_single = Y_single - Y_ref;
err_fixed  = Y_fixed_dbl - Y_ref;

fprintf('---------------------------------------------\n');
fprintf('FFT Vector Length        : %d\n', Npts);
fprintf('Q Format                 : Q%d.%d\n', x, y);
fprintf('Max Abs Error (double)   : %e\n', max(abs(err_double)));
fprintf('Mean Abs Error (double)  : %e\n', mean(abs(err_double)));
fprintf('Max Abs Error (single)   : %e\n', max(abs(err_single)));
fprintf('Mean Abs Error (single)  : %e\n', mean(abs(err_single)));
fprintf('Max Abs Error (fixed)    : %e\n', max(abs(err_fixed)));
fprintf('Mean Abs Error (fixed)   : %e\n', mean(abs(err_fixed)));
fprintf('---------------------------------------------\n');

disp('Index |  A (real + j imag) |  B (real + j imag) |      HW Double     |      HW Single     |      HW Fixed      |        Reference');
disp('--------------------------------------------------------------------------------------------------------------------------------------');

for k = 1:5
    fprintf('%4d  | %+1.5f %+1.5fj | %+1.5f %+1.5fj | %+1.5f %+1.5fj | %+1.5f %+1.5fj | %+1.5f %+1.5fj | %+1.5f %+1.5fj\n', ...
        k, ...
        real(A(k)), imag(A(k)), ...
        real(B(k)), imag(B(k)), ...
        real(Y_double(k)), imag(Y_double(k)), ...
        real(Y_single(k)), imag(Y_single(k)), ...
        real(Y_fixed_dbl(k)), imag(Y_fixed_dbl(k)), ...
        real(Y_ref(k)), imag(Y_ref(k)));
end

