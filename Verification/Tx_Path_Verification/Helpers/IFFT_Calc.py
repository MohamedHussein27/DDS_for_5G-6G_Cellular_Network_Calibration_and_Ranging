import numpy as np
import matplotlib.pyplot as plt

def simulate_and_visualize_ifft(n_fft=64, mode="directed", directed_real=None, directed_imag=None):
    """
    Calculates, prints, and plots the IFFT.
    
    Modes:
      - "directed": Uses exact user-provided real and imaginary vectors (e.g., FFT outputs).
      - "random": Generates continuous complex values to mimic raw FFT/OFDM outputs.
    """
    
    if mode == "directed":
        if directed_real is None or directed_imag is None:
            raise ValueError("Directed mode requires both real and imag vectors.")
        if len(directed_real) != n_fft or len(directed_imag) != n_fft:
            raise ValueError(f"Vectors must match n_fft size ({n_fft}).")
            
        freq_symbols = np.array(directed_real, dtype=float) + 1j * np.array(directed_imag, dtype=float)
        print(f"=== Running DIRECTED Test (N={n_fft}) ===")

    elif mode == "random":
        # Generate complex Gaussian data to mimic the variable magnitudes of an FFT output
        np.random.seed(42)
        real_parts = np.random.normal(0, 1, size=n_fft)
        imag_parts = np.random.normal(0, 1, size=n_fft)
        freq_symbols = real_parts + 1j * imag_parts
        
        # Apply DC null (Bin 0)
        freq_symbols[0] = 0 
        print(f"=== Running RANDOM (Variable Magnitude) Test (N={n_fft}) ===")
        
    else:
        raise ValueError("Invalid mode.")

    # --- 1. CALCULATE IFFT ---
    time_samples = np.fft.ifft(freq_symbols)
    
    # --- 2. CALCULATE AND PRINT METRICS ---
    energy_freq = np.sum(np.abs(freq_symbols)**2) / n_fft
    energy_time = np.sum(np.abs(time_samples)**2)
    
    power_time = np.abs(time_samples)**2
    peak_power = np.max(power_time)
    avg_power = np.mean(power_time)
    papr_db = 10 * np.log10(peak_power / avg_power)
    
    print(f"Energy (Freq / N):  {energy_freq:.6f}")
    print(f"Energy (Time):      {energy_time:.6f}")
    print(f"Energy Match:       {np.isclose(energy_freq, energy_time)}")
    print(f"PAPR:               {papr_db:.2f} dB\n")
    
    print("First 8 Time-Domain Samples:")
    for i in range(min(8, n_fft)):
        print(f"[{i}]: {time_samples[i].real:+.6f} {time_samples[i].imag:+.6f}j")

    # --- 3. VISUALIZATION ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot 1: Frequency Domain Input (Notice the variable magnitudes now)
    ax1.stem(range(n_fft), np.abs(freq_symbols), basefmt="b-")
    ax1.set_title("Frequency Domain Input (Magnitude of FFT Output / OFDM Data)")
    ax1.set_xlabel("Frequency Bin")
    ax1.set_ylabel("Magnitude")
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Plot 2: Time Domain Output
    ax2.plot(range(n_fft), time_samples.real, label="Real (I)", marker='.', markersize=6)
    ax2.plot(range(n_fft), time_samples.imag, label="Imaginary (Q)", marker='x', markersize=4, alpha=0.7)
    ax2.set_title("Time Domain Output (IFFT Result)")
    ax2.set_xlabel("Sample Index")
    ax2.set_ylabel("Amplitude")
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.show()

    return time_samples

# ==========================================
# EXECUTION EXAMPLE: Feeding Variable Data
# ==========================================
FFT_SIZE = 4096

# Example: We simulate an upstream FFT output where subcarriers have varying amplitudes
# (e.g., due to channel fading or radar target reflections)
#custom_real = np.random.uniform(-5.5, 5.5, FFT_SIZE)
#custom_imag = np.random.uniform(-5.5, 5.5, FFT_SIZE)

custom_real = [ 0,  1, -1,  1,  0,  0, -1,  1]
custom_imag = [ 0,  1,  1, -1,  0,  0, -1,  0]


# Directed test with a single impulse at bin 1 (index 1)
custom_real = [0, 1] + [0] * (FFT_SIZE - 2)   # index0=0, index1=1, rest zeros
custom_imag = [0] * FFT_SIZE                  # all zeros


# Manually insert your specific Guardband (e.g., bins mapped to 200kHz-210kHz)
# For this example, we assume bins 12 to 15 represent that guardband
#guardband_start = 12
#guardband_end = 15
#custom_real[guardband_start:guardband_end+1] = 0
#custom_imag[guardband_start:guardband_end+1] = 0

out_samples = simulate_and_visualize_ifft(
    n_fft=FFT_SIZE, 
    mode="directed", 
    directed_real=custom_real, 
    directed_imag=custom_imag
)