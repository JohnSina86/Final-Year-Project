import numpy as np
import cupy as cp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import time
import sys
import os

# Accept directory argument
if len(sys.argv) > 1:
    os.chdir(sys.argv[1])

# Physical constants
MHz = np.float32(1e6)
f = np.float32(2000.0) * MHz
c0 = np.float32(3e8)
omega = np.float32(2) * np.pi * f
k0 = omega / c0
lambda_val = c0 / f
eta = np.float32(377.0)

print(f"Frequency: {f/MHz:.0f} MHz")
print(f"Wavelength: {lambda_val:.4f} m")
print(f"Wavenumber k0: {k0:.2f} rad/m")
print(f"GPU: {cp.cuda.Device().name}")

def read_geometry(filename='geometry.txt'):
    """Read segment positions"""
    positions = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('%') or line.strip() == '':
                continue
            x, y = map(float, line.split())
            positions.append((x, y))
    return np.array(positions, dtype=np.float32)

def hankel2_approx_gpu(z):
    """GPU-accelerated Hankel function H0(2)"""
    result = cp.zeros(z.shape, dtype=cp.complex64)
    
    small_mask = z < 0.01
    large_mask = ~small_mask
    
    if cp.any(small_mask):
        z_small = z[small_mask]
        result[small_mask] = 1.0 - 1j * (2/cp.pi) * cp.log(z_small/2 + 1e-10)
    
    if cp.any(large_mask):
        z_large = z[large_mask]
        sqrt_term = cp.sqrt(2 / (cp.pi * z_large))
        exp_term = cp.exp(-1j * (z_large + cp.pi / 4))
        result[large_mask] = sqrt_term * exp_term
    
    return result

def main():
    start_total = time.time()
    
    # Read parameters
    try:
        with open('Transmitter pos', 'r') as f:
            strip_length = float(f.readline())
            TX = np.float32(float(f.readline()))
            TY = np.float32(float(f.readline()))
        print(f"Strip length: {strip_length:.3f} m")
        print(f"Transmitter: ({TX:.2f}, {TY:.2f})")
    except FileNotFoundError:
        print("ERROR: 'Transmitter pos' file not found!")
        return
    
    # Read geometry
    try:
        positions = read_geometry()
        x_positions = positions[:, 0]
        y_positions = positions[:, 1]
        print(f"Loaded {len(positions)} segment positions")
    except FileNotFoundError:
        print("ERROR: 'geometry.txt' file not found!")
        return
    
    # Read current coefficients
    try:
        x_current_real = []
        x_current_imag = []
        with open('current_distribution.txt', 'r') as f:
            for line in f:
                if line.startswith('%') or not line.strip():
                    continue
                parts = line.strip().replace('(', '').replace(')', '').split(',')
                if len(parts) == 2:
                    x_current_real.append(float(parts[0]))
                    x_current_imag.append(float(parts[1]))
        
        x_current = (np.array(x_current_real, dtype=np.float32) + 
                     1j * np.array(x_current_imag, dtype=np.float32)).astype(np.complex64)
        N = len(x_current)
        print(f"Loaded {N} current coefficients")
        
        if N == 0:
            print("ERROR: No current data!")
            return
    except FileNotFoundError:
        print("ERROR: 'current_distribution.txt' not found!")
        return
    
    # Grid size
    grid_size = 1500
    print(f"\nGrid size: {grid_size}x{grid_size}")
    
    # Discretization
    delta_s = np.float32(lambda_val / 40)
    
    # Observation grid
    x_obs = np.linspace(-0.1, 0.6, grid_size, dtype=np.float32)
    y_obs = np.linspace(-0.1, 0.6, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x_obs, y_obs, copy=False)
    
    print(f"Computing fields on {grid_size}x{grid_size} grid...")
    
    # Transfer to GPU
    X_gpu = cp.asarray(X)
    Y_gpu = cp.asarray(Y)
    
    eta_k0_4 = eta * k0 / 4
    
    # ===== INCIDENT FIELD (GPU) =====
    print("\n[GPU] Computing incident field...")
    start_inc = time.time()
    
    dx_inc = X_gpu - TX
    dy_inc = Y_gpu - TY
    R_inc = cp.sqrt(dx_inc**2 + dy_inc**2)
    R_inc = cp.where(R_inc < 1e-6, 1e-6, R_inc)
    
    E_field_gpu = hankel2_approx_gpu(k0 * R_inc)
    
    del dx_inc, dy_inc, R_inc
    
    print(f"  Time: {time.time() - start_inc:.2f}s")
    
    # ===== SCATTERED FIELD (GPU) =====
    print(f"\n[GPU] Computing scattered field...")
    start_scatter = time.time()
    
    # Transfer source data to GPU
    x_src_gpu = cp.asarray(x_positions)
    y_src_gpu = cp.asarray(y_positions)
    I_gpu = cp.asarray(x_current)
    
    # Process in batches
    batch_size = 50
    
    for batch_start in range(0, N, batch_size):
        batch_end = min(batch_start + batch_size, N)
        
        for i in range(batch_start, batch_end):
            dx = X_gpu - x_src_gpu[i]
            dy = Y_gpu - y_src_gpu[i]
            R = cp.sqrt(dx**2 + dy**2)
            R = cp.where(R < delta_s / 10, delta_s / 10, R)
            
            G = hankel2_approx_gpu(k0 * R)
            E_field_gpu += G * I_gpu[i]
        
        # Progress
        if (batch_end % 100) == 0 or batch_end == N:
            progress = batch_end / N * 100
            elapsed = time.time() - start_scatter
            eta_time = elapsed / batch_end * (N - batch_end) if batch_end > 0 else 0
            print(f"  Progress: {progress:.0f}% ({batch_end}/{N}), ETA: {eta_time:.1f}s")
    
    # Apply factor
    E_field_gpu *= -(eta_k0_4 * delta_s)
    
    print(f"  GPU time: {time.time() - start_scatter:.2f}s")
    
    # Transfer back to CPU
    print("Transferring to CPU...")
    E_field = cp.asnumpy(E_field_gpu)
    E_mag = np.abs(E_field)
    
    total_time = time.time() - start_total
    print(f"\nField range: [{E_mag.min():.4e}, {E_mag.max():.4e}] V/m")
    print(f"Total time: {total_time:.2f}s")
    
    # Plotting
    print("Generating plot...")
    colors = ['#FFFFFF', '#E0FFFF', '#B0E0E6', '#87CEEB', '#4169E1', '#0000CD', '#00008B']
    cmap_custom = LinearSegmentedColormap.from_list('white_to_blue', colors, N=256)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    extent = [x_obs[0], x_obs[-1], y_obs[0], y_obs[-1]]
    im = ax.imshow(E_mag, extent=extent, origin='lower', cmap=cmap_custom, 
                   aspect='equal', interpolation='bilinear')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.invert_yaxis()
    cbar.set_ticks([E_mag.min(), E_mag.max()])
    cbar.set_ticklabels(['Strongest', 'Weakest'])
    
    downsample = max(1, len(x_positions) // 200)
    ax.plot(x_positions[::downsample], y_positions[::downsample], 
            'c-', linewidth=2, label='L-bracket', alpha=0.8)
    ax.scatter(x_positions[::downsample], y_positions[::downsample], 
               c='cyan', s=30, edgecolors='white', linewidth=0.5, zorder=10, alpha=0.9)
    
    ax.set_xlabel('X Position (m)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y Position (m)', fontsize=12, fontweight='bold')
    ax.set_title(f'GPU Field ({grid_size}×{grid_size})', fontsize=14, fontweight='bold')
    ax.set_xlim(-0.1, 0.6)
    ax.set_ylim(-0.1, 0.6)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, color='gray', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig('total_fields_gpu.png', dpi=200, bbox_inches='tight')
    print(f"Saved: total_fields_gpu.png")
    print(f"🚀 Done! {total_time:.1f}s")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
