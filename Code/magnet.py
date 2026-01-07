import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from multiprocessing import Pool, cpu_count
import time

# Physical constants (float32)
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

def hankel2_approx_vec(z):
    """Vectorized Hankel function"""
    result = np.zeros_like(z, dtype=np.complex64)
    small_mask = np.abs(z) < 0.01
    large_mask = ~small_mask
    
    if np.any(small_mask):
        result[small_mask] = 1.0 - 1j * (2/np.pi) * np.log(z[small_mask]/2)
    
    if np.any(large_mask):
        z_large = z[large_mask]
        result[large_mask] = np.sqrt(2 / (np.pi * z_large)) * np.exp(-1j * (z_large + np.pi / 4))
    
    return result

def compute_scattered_chunk(args):
    """Compute scattered field for a grid chunk - multiprocessing worker"""
    X_chunk, Y_chunk, x_src, y_src, x_current, delta_s, eta_k0_4, k0, chunk_id = args
    
    chunk_shape = X_chunk.shape
    N = len(x_src)
    
    # Reshape for broadcasting
    X_flat = X_chunk[:, :, np.newaxis]
    Y_flat = Y_chunk[:, :, np.newaxis]
    
    x_src_bc = x_src[np.newaxis, np.newaxis, :]
    y_src_bc = y_src[np.newaxis, np.newaxis, :]
    
    # Distance calculation
    dx = X_flat - x_src_bc
    dy = Y_flat - y_src_bc
    R = np.sqrt(dx**2 + dy**2)
    
    # Avoid singularities
    R[R < delta_s / 10] = delta_s / 10
    
    # Green's function
    G = hankel2_approx_vec(k0 * R)
    
    # Current contribution
    I = x_current[np.newaxis, np.newaxis, :]
    
    # Scattered field
    E_scattered = -(eta_k0_4 * delta_s) * np.sum(G * I, axis=2)
    
    return chunk_id, E_scattered

def main():
    start_total = time.time()
    
    # Read parameters
    with open('Transmitter pos', 'r') as f:
        strip_length = float(f.readline())
        TX = np.float32(float(f.readline()))
        TY = np.float32(float(f.readline()))
    
    print(f"Strip length: {strip_length:.3f} m")
    print(f"Transmitter: ({TX:.2f}, {TY:.2f})")
    
    # Read current coefficients (float32)
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
    
    x_current = np.array(x_current_real, dtype=np.float32) + 1j * np.array(x_current_imag, dtype=np.float32)
    N = len(x_current)
    print(f"Loaded {N} current coefficients from current_distribution.txt")
    
    # Strip positions
    delta_s = np.float32(lambda_val / 10)
    x_positions = np.array([(i + 0.5) * delta_s for i in range(N)], dtype=np.float32)
    y_positions = np.zeros(N, dtype=np.float32)
    
    print(f"Strip: x=[{x_positions[0]:.3f}, {x_positions[-1]:.3f}] m, y=0")
    
    # Observation grid (float32)
    grid_size = 2500
    x_obs = np.linspace(-2.0, 2.5, grid_size, dtype=np.float32)
    y_obs = np.linspace(-2.0, 2.0, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x_obs, y_obs)
    
    print(f"\nComputing fields on {grid_size}x{grid_size} grid...")
    
    # Constants
    eta_k0_4 = eta * k0 / 4
    
    # ===== INCIDENT FIELD =====
    print("Computing incident field...")
    start_inc = time.time()
    dx_inc = X - TX
    dy_inc = Y - TY
    R_inc = np.sqrt(dx_inc**2 + dy_inc**2)
    R_inc[R_inc < 1e-6] = 1e-6
    
    E_field = hankel2_approx_vec(k0 * R_inc).astype(np.complex64)
    print(f"  Incident field: {time.time() - start_inc:.2f}s")
    
    # ===== MULTIPROCESSING SCATTERED FIELD =====
    print(f"Computing scattered field (multiprocessing, {cpu_count()} cores)...")
    start_scatter = time.time()
    
    # Chunk parameters
    chunk_size = 250  # Smaller chunks for better load balancing
    num_chunks_y = (grid_size + chunk_size - 1) // chunk_size
    num_chunks_x = (grid_size + chunk_size - 1) // chunk_size
    
    # Prepare chunks
    chunks = []
    chunk_positions = []
    
    for chunk_i in range(num_chunks_y):
        for chunk_j in range(num_chunks_x):
            row_start = chunk_i * chunk_size
            row_end = min(row_start + chunk_size, grid_size)
            col_start = chunk_j * chunk_size
            col_end = min(col_start + chunk_size, grid_size)
            
            X_chunk = X[row_start:row_end, col_start:col_end]
            Y_chunk = Y[row_start:row_end, col_start:col_end]
            
            chunk_id = (row_start, row_end, col_start, col_end)
            chunks.append((X_chunk, Y_chunk, x_positions, y_positions, 
                          x_current, delta_s, eta_k0_4, k0, chunk_id))
            chunk_positions.append(chunk_id)
    
    # Process chunks in parallel
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(compute_scattered_chunk, chunks)
    
    # Accumulate results
    for chunk_id, E_scatter in results:
        row_start, row_end, col_start, col_end = chunk_id
        E_field[row_start:row_end, col_start:col_end] += E_scatter
    
    print(f"  Scattered field: {time.time() - start_scatter:.2f}s")
    
    # Magnitude
    E_mag = np.abs(E_field)
    print(f"\nField range: [{E_mag.min():.4e}, {E_mag.max():.4e}] V/m")
    print(f"Total computation time: {time.time() - start_total:.2f}s")
    
    # Plotting
    colors = ['#FFFFFF', '#E0FFFF', '#B0E0E6', '#87CEEB', '#4169E1', '#0000CD', '#00008B']
    n_bins = 1048576
    cmap_custom = LinearSegmentedColormap.from_list('white_to_blue', colors, N=n_bins)
    
    fig, ax = plt.subplots(figsize=(8, 10))
    im = ax.pcolormesh(X, Y, E_mag, shading='auto', cmap=cmap_custom)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.invert_yaxis()
    cbar.set_ticks([E_mag.min(), E_mag.max()])
    cbar.set_ticklabels(['Strongest', 'Weakest'])
    
    ax.plot(x_positions, y_positions, 'k-', linewidth=3, label='Metal strip')
    ax.set_xlabel('x position (m)', fontsize=12)
    ax.set_ylabel('y position (m)', fontsize=12)
    ax.set_title('Total Electric Field Magnitude (Multiprocessing)', fontsize=14)
    ax.set_xlim(-2, 2.5)
    ax.set_ylim(-2, 2)
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, color='gray', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig('total_fields.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved as total_fields.png")
    plt.show()

if __name__ == "__main__":
    main()
