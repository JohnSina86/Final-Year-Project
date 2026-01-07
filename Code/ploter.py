import numpy as np
import matplotlib.pyplot as plt
import os

# Set working directory to script location
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

print(f"Working directory: {os.getcwd()}")

# Read solution from new file
solution = []
with open('current_distribution.txt', 'r') as f:
    for line in f:
        line = line.strip()

        # Skip comment lines
        if line.startswith('%') or line == '':
            continue

        # Parse (real,imag)
        real, imag = line.strip('()').split(',')
        solution.append(complex(float(real), float(imag)))

solution = np.array(solution)
print(f"Read {len(solution)} values from current_distribution.txt")


# Generate positions based on problem size
N = len(solution)
# Assuming strip starts at 0 and uses lambda/10 discretization
# For 2 GHz, lambda = c/f = 3e8/2e9 = 0.15 m
# With length_of_strip = 3.0, delta_s = 0.15/10 = 0.015
# positions are (i+0.5)*delta_s for i = 0 to N-1

# Calculate from N
lambda_val = 0.15  # wavelength at 2 GHz in meters
delta_s = lambda_val / 10  # discretization
positions = np.array([(i + 0.5) * delta_s for i in range(N)])

print(f"Generated {len(positions)} positions")

# Calculate magnitude
mag_solution = np.abs(solution)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(positions, mag_solution, 'b-', linewidth=2, label='Current Distribution (SSOR)')

plt.xlabel('distance on plate', fontsize=12)
plt.ylabel('abs(current)', fontsize=12)
plt.title('Current on plate after SSOR iterations', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save the figure
plt.savefig('current_distribution_plot.png', dpi=300, bbox_inches='tight')
print("Plot saved as current_distribution_plot.png")

plt.show()
