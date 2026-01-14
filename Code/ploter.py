import numpy as np
import matplotlib.pyplot as plt
import os
import re


# Set working directory to script location
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)


print(f"Working directory: {os.getcwd()}")


def read_geometry(filename='geometry.txt'):
    """Read segment positions from geometry file"""
    positions = []
    N_strip1 = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('%'):
                if 'N_horizontal' in line:
                    match = re.search(r'N_horizontal[=\s]+(\d+)', line)
                    if match:
                        N_strip1 = int(match.group(1))
                        print(f"Parsed N_horizontal = {N_strip1}")
                continue
            if line == '':
                continue

            parts = line.split()
            if len(parts) >= 2:
                x, y = float(parts[0]), float(parts[1])
                positions.append((x, y))

    # Fallback estimation
    if N_strip1 is None:
        N_strip1 = sum(1 for p in positions if abs(p[1]) < 1e-6)
        print(f"Estimated N_horizontal = {N_strip1}")

    return positions, N_strip1


def read_currents(filename='current_distribution.txt'):
    """Read current distribution from file"""
    currents = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('%') or line == '':
                continue
            line_clean = line.strip('()')
            parts = line_clean.split(',')
            if len(parts) == 2:
                real = float(parts[0])
                imag = float(parts[1])
                currents.append(complex(real, imag))
    return np.array(currents)


# Read geometry and currents
print("\nReading files...")
positions, N_strip1 = read_geometry()
currents = read_currents()

N_total = len(currents)
N_strip2 = N_total - N_strip1

print(f"Total segments: {N_total} (Horizontal: {N_strip1}, Vertical: {N_strip2})")

# Extract positions
x_pos = np.array([p[0] for p in positions])
y_pos = np.array([p[1] for p in positions])

# Current magnitudes and phases
current_mag = np.abs(currents)
current_phase = np.angle(currents)

# Create figure with 2 subplots side-by-side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5), dpi=100)

# =================================================================
# Plot 1: Current magnitude vs segment index
# =================================================================
ax1.plot(range(N_strip1), current_mag[:N_strip1], 
         'r-', linewidth=2, alpha=0.8, label='Horizontal arm')
ax1.plot(range(N_strip1, N_total), current_mag[N_strip1:], 
         'b-', linewidth=2, alpha=0.8, label='Vertical arm')
ax1.axvline(N_strip1-0.5, color='gray', linestyle='--', linewidth=1.5, 
            alpha=0.5, label='Junction')
ax1.set_xlabel('Segment Index', fontsize=12, fontweight='bold')
ax1.set_ylabel('Current Magnitude (A)', fontsize=12, fontweight='bold')
ax1.set_title('Current Magnitude Distribution', fontsize=14, fontweight='bold')
ax1.legend(fontsize=11, framealpha=0.9)
ax1.grid(True, alpha=0.3, linestyle=':', linewidth=1)
ax1.tick_params(labelsize=11)

# =================================================================
# Plot 2: Current phase vs segment index
# =================================================================
ax2.plot(range(N_strip1), current_phase[:N_strip1], 
         'r-', linewidth=2, alpha=0.8, label='Horizontal arm')
ax2.plot(range(N_strip1, N_total), current_phase[N_strip1:], 
         'b-', linewidth=2, alpha=0.8, label='Vertical arm')
ax2.axvline(N_strip1-0.5, color='gray', linestyle='--', linewidth=1.5, 
            alpha=0.5, label='Junction')
ax2.set_xlabel('Segment Index', fontsize=12, fontweight='bold')
ax2.set_ylabel('Phase (radians)', fontsize=12, fontweight='bold')
ax2.set_title('Current Phase Distribution', fontsize=14, fontweight='bold')
ax2.legend(fontsize=11, framealpha=0.9)
ax2.grid(True, alpha=0.3, linestyle=':', linewidth=1)
ax2.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig('current_distribution_plot.png', dpi=150, bbox_inches='tight')
print(f"\nPlot saved: current_distribution_plot.png")
print(f"Max current magnitude: {current_mag.max()*1e3:.4f} mA")
print(f"Min current magnitude: {current_mag.min()*1e3:.4f} mA")
plt.show()
