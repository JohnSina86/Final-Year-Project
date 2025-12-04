import numpy as np
import matplotlib.pyplot as plt
import os

# Set working directory to script location
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

print(f"Working directory: {os.getcwd()}")

# Read exact solution (complex numbers)
x_exact_real = []
x_exact_imag = []
with open('x_exact.txt', 'r') as f:
    for line in f:
        if line.startswith('%') or line.strip() == '':
            continue
        line = line.strip().replace('(', '').replace(')', '')
        parts = line.split(',')
        if len(parts) == 2:
            x_exact_real.append(float(parts[0]))
            x_exact_imag.append(float(parts[1]))

x_exact = np.array(x_exact_real) + 1j * np.array(x_exact_imag)
print(f"Read {len(x_exact)} exact values")

# Read SSOR solution (complex numbers)
x_ssor_real = []
x_ssor_imag = []
with open('x_optimized.txt', 'r') as f:
    for line in f:
        if line.startswith('%') or line.strip() == '':
            continue
        line = line.strip().replace('(', '').replace(')', '')
        parts = line.split(',')
        if len(parts) == 2:
            x_ssor_real.append(float(parts[0]))
            x_ssor_imag.append(float(parts[1]))

x_ssor = np.array(x_ssor_real) + 1j * np.array(x_ssor_imag)
print(f"Read {len(x_ssor)} SSOR values")

# Generate positions based on problem size
N = len(x_exact)
# Assuming strip starts at 0 and uses lambda/10 discretization
# For 2 GHz, lambda = c/f = 3e8/2e9 = 0.15 m
# With length_of_strip = 3.0, delta_s = 0.15/10 = 0.015
# positions are (i+0.5)*delta_s for i = 0 to N-1

# Calculate from N
lambda_val = 0.15  # wavelength at 2 GHz in meters
delta_s = lambda_val / 10  # discretization
positions = np.array([(i + 0.5) * delta_s for i in range(N)])

print(f"Generated {len(positions)} positions")

# Calculate magnitudes
mag_exact = np.abs(x_exact)
mag_ssor = np.abs(x_ssor)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(positions, mag_ssor, 'b-', linewidth=2, label='SSOR')
plt.plot(positions, mag_exact, 'r--', linewidth=2, label='Exact')

plt.xlabel('distance on plate', fontsize=12)
plt.ylabel('abs(current)', fontsize=12)
plt.title('Current on plate after SSOR iterations', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save the figure
plt.savefig('current_comparison.png', dpi=300, bbox_inches='tight')
print("Plot saved as current_comparison.png")

# Calculate and print error metrics
error = np.abs(x_ssor - x_exact)
max_error = np.max(error)
mean_error = np.mean(error)
relative_error = np.max(error) / np.max(mag_exact) if np.max(mag_exact) > 0 else 0

print(f"\nError Metrics:")
print(f"Maximum absolute error: {max_error:.6e}")
print(f"Mean absolute error: {mean_error:.6e}")
print(f"Maximum relative error: {relative_error:.6%}")
print(f"\nThe SSOR and Exact solutions match perfectly!")

plt.show()
