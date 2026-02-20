import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from multiprocessing import Pool, cpu_count
import time
import sys
import psutil
import os

# Physical constants (float32 for memory efficiency)
MHz = np.float32(1e6)
try:
    with open("frequency.txt", "r") as _ff:
        f = np.float32(float(_ff.read().strip()) * 1e9)
except FileNotFoundError:
    f = np.float32(2000.0 * 1e6)
c0 = np.float32(3e8)
omega = np.float32(2) * np.pi * f
k0 = omega / c0
lambda_val = c0 / f
eta = np.float32(377.0)

print(f"Frequency: {f/MHz:.0f} MHz")
print(f"Wavelength: {lambda_val:.4f} m")
print(f"Wavenumber k0: {k0:.2f} rad/m")

def get_memory_usage():
    """Get current memory usage in MB"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

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

def hankel2_approx_vec(z):
    """Optimized vectorized Hankel function with memory management"""
    z = np.asarray(z, dtype=np.float32)
    result = np.zeros(z.shape, dtype=np.complex64)

    small_mask = z < 0.01
    large_mask = ~small_mask

    # Process in chunks to avoid memory issues
    if np.any(small_mask):
        z_small = z[small_mask]
        result[small_mask] = 1.0 - 1j * (2/np.pi) * np.log(z_small/2 + 1e-10)

    if np.any(large_mask):
        z_large = z[large_mask]
        sqrt_term = np.sqrt(2 / (np.pi * z_large))
        exp_term = np.exp(-1j * (z_large + np.pi / 4))
        result[large_mask] = sqrt_term * exp_term

    return result

def compute_scattered_chunk_optimized(args):
    """OPTIMIZED: Process smaller segments to avoid memory overflow"""
    X_chunk, Y_chunk, x_src, y_src, x_current, delta_s, eta_k0_4, k0_val, chunk_id = args

    ny, nx = X_chunk.shape
    N_src = len(x_src)
    E_scattered = np.zeros((ny, nx), dtype=np.complex64)

    # Process sources in smaller batches to control memory
    batch_size = min(50, N_src)  # Process 50 sources at a time

    for batch_start in range(0, N_src, batch_size):
        batch_end = min(batch_start + batch_size, N_src)

        # Extract batch
        x_batch = x_src[batch_start:batch_end]
        y_batch = y_src[batch_start:batch_end]
        I_batch = x_current[batch_start:batch_end]

        # Compute distances for this batch
        for i in range(len(x_batch)):
            dx = X_chunk - x_batch[i]
            dy = Y_chunk - y_batch[i]
            R = np.sqrt(dx**2 + dy**2)
            R[R < delta_s / 10] = delta_s / 10

            # Green's function
            G = hankel2_approx_vec(k0_val * R)

            # Accumulate contribution
            E_scattered += G * I_batch[i]

    # Apply factor
    E_scattered *= -(eta_k0_4 * delta_s)

    return chunk_id, E_scattered

def adaptive_grid_size(N_segments, available_memory_gb):
    """Determine safe grid size based on data size and available memory"""
    # Estimate memory per grid point (complex64 + processing overhead)
    bytes_per_point = 32  # Conservative estimate

    # Reserve memory for arrays and processing
    available_bytes = available_memory_gb * 1e9 * 0.5  # Use 50% of available

    # Calculate maximum grid points
    max_points = int(available_bytes / bytes_per_point)
    max_grid_size = int(np.sqrt(max_points))

    # Cap between reasonable limits
    min_size = 1750
    max_size = 3000

    # Adjust based on number of segments
    if N_segments < 100:
        suggested = 3000
    elif N_segments < 300:
        suggested = 2500
    elif N_segments < 600:
        suggested = 2000
    else:
        suggested = 1750

    grid_size = min(max(suggested, min_size), max_size, max_grid_size)

    print(f"  Segments: {N_segments}")
    print(f"  Available memory: {available_memory_gb:.1f} GB")
    print(f"  Suggested grid size: {grid_size}x{grid_size}")

    return grid_size

def main():
    start_total = time.time()

    print(f"\nInitial memory usage: {get_memory_usage():.1f} MB")

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
            print("ERROR: No current data loaded!")
            return

    except FileNotFoundError:
        print("ERROR: 'current_distribution.txt' file not found!")
        return

    # Check available memory
    available_memory = psutil.virtual_memory().available / (1024**3)  # GB
    total_memory = psutil.virtual_memory().total / (1024**3)
    print(f"\nMemory: {available_memory:.1f} GB available / {total_memory:.1f} GB total")

    # Adaptive grid size
    grid_size = adaptive_grid_size(N, available_memory)

    # Discretization
    delta_s = np.float32(lambda_val / 100)  # Must match SEGMENTS_PER_LAMBDA in mom.c

    # Observation grid
    x_obs = np.linspace(-0.1, 0.6, grid_size, dtype=np.float32)
    y_obs = np.linspace(-0.1, 0.6, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x_obs, y_obs, copy=False)  # Don't copy unnecessarily

    print(f"\nComputing fields on {grid_size}x{grid_size} grid...")
    print(f"Memory after grid creation: {get_memory_usage():.1f} MB")

    # Constants
    eta_k0_4 = eta * k0 / 4

    # ===== INCIDENT FIELD =====
    print("\nComputing incident field...")
    start_inc = time.time()

    dx_inc = X - TX
    dy_inc = Y - TY
    R_inc = np.sqrt(dx_inc**2 + dy_inc**2, dtype=np.float32)
    R_inc[R_inc < 1e-6] = 1e-6

    E_field = hankel2_approx_vec(k0 * R_inc)

    # Free memory
    del dx_inc, dy_inc, R_inc

    print(f"  Time: {time.time() - start_inc:.2f}s")
    print(f"  Memory: {get_memory_usage():.1f} MB")

    # ===== MULTIPROCESSING SCATTERED FIELD =====
    n_cores = max(1, cpu_count() - 1)  # Leave one core free
    print(f"\nComputing scattered field ({n_cores} cores)...")
    start_scatter = time.time()

    # Adaptive chunk size based on grid size
    if grid_size <= 800:
        chunk_size = 200
    elif grid_size <= 1500:
        chunk_size = 250
    else:
        chunk_size = 300

    num_chunks_y = (grid_size + chunk_size - 1) // chunk_size
    num_chunks_x = (grid_size + chunk_size - 1) // chunk_size
    total_chunks = num_chunks_y * num_chunks_x

    print(f"  Chunk size: {chunk_size}x{chunk_size}")
    print(f"  Total chunks: {total_chunks}")

    # Prepare chunks
    chunks = []
    for chunk_i in range(num_chunks_y):
        for chunk_j in range(num_chunks_x):
            row_start = chunk_i * chunk_size
            row_end = min(row_start + chunk_size, grid_size)
            col_start = chunk_j * chunk_size
            col_end = min(col_start + chunk_size, grid_size)

            # Extract chunk views (no copy)
            X_chunk = X[row_start:row_end, col_start:col_end].copy()  # Copy for multiprocessing
            Y_chunk = Y[row_start:row_end, col_start:col_end].copy()

            chunk_id = (row_start, row_end, col_start, col_end)
            chunks.append((X_chunk, Y_chunk, x_positions, y_positions, 
                          x_current, delta_s, eta_k0_4, k0, chunk_id))

    print(f"  Processing {total_chunks} chunks...")

    # Process chunks in parallel with progress
    try:
        with Pool(processes=n_cores) as pool:
            results = []
            for i, result in enumerate(pool.imap(compute_scattered_chunk_optimized, chunks)):
                results.append(result)
                if (i + 1) % max(1, total_chunks // 10) == 0:
                    progress = (i + 1) / total_chunks * 100
                    elapsed = time.time() - start_scatter
                    eta_time = elapsed / (i + 1) * (total_chunks - i - 1)
                    print(f"    Progress: {progress:.0f}% ({i+1}/{total_chunks}), "
                          f"ETA: {eta_time:.1f}s")
    except Exception as e:
        print(f"ERROR during multiprocessing: {e}")
        print("Trying single-threaded fallback...")
        results = []
        for i, chunk in enumerate(chunks):
            results.append(compute_scattered_chunk_optimized(chunk))
            if (i + 1) % max(1, total_chunks // 10) == 0:
                print(f"    Progress: {(i+1)/total_chunks*100:.0f}% ({i+1}/{total_chunks})")

    # Accumulate results
    print("  Assembling results...")
    for chunk_id, E_scatter in results:
        row_start, row_end, col_start, col_end = chunk_id
        E_field[row_start:row_end, col_start:col_end] += E_scatter

    print(f"  Time: {time.time() - start_scatter:.2f}s")
    print(f"  Memory: {get_memory_usage():.1f} MB")

    # Magnitude
    E_mag = np.abs(E_field)
    print(f"\nField range: [{E_mag.min():.4e}, {E_mag.max():.4e}] V/m")
    print(f"Total computation time: {time.time() - start_total:.2f}s")
    print(f"Final memory usage: {get_memory_usage():.1f} MB")

    # Plotting
    print("\nGenerating plot...")
    colors = ['#FFFFFF', '#E0FFFF', '#B0E0E6', '#87CEEB', '#4169E1', '#0000CD', '#00008B']
    n_bins = 256  # Reduced from 1048576 - no visual difference
    cmap_custom = LinearSegmentedColormap.from_list('white_to_blue', colors, N=n_bins)

    fig, ax = plt.subplots(figsize=(10, 10))

    # Use imshow for better performance on large grids
    extent = [x_obs[0], x_obs[-1], y_obs[0], y_obs[-1]]
    im = ax.imshow(E_mag, extent=extent, origin='lower', cmap=cmap_custom, 
                   aspect='equal', interpolation='bilinear')

    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.invert_yaxis()
    cbar.set_ticks([E_mag.min(), E_mag.max()])
    cbar.set_ticklabels(['Strongest', 'Weakest'])

    # Downsample L-bracket visualization for large datasets
    downsample = max(1, len(x_positions) // 200)
    ax.plot(x_positions[::downsample], y_positions[::downsample], 
            'c-', linewidth=2, label='L-bracket', alpha=0.8)
    ax.scatter(x_positions[::downsample], y_positions[::downsample], 
               c='cyan', s=30, edgecolors='white', linewidth=0.5, zorder=10, alpha=0.9)

    ax.set_xlabel('X Position (m)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y Position (m)', fontsize=12, fontweight='bold')
    ax.set_title(f'Total Electric Field Magnitude ({grid_size}×{grid_size} grid)', 
                 fontsize=14, fontweight='bold')
    ax.set_xlim(-0.1, 0.6)
    ax.set_ylim(-0.1, 0.6)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, color='gray', linewidth=0.5)

    plt.tight_layout()

    # Adaptive DPI based on grid size
    dpi = min(300, max(150, 30000 // grid_size))
    plt.savefig('total_fields.png', dpi=dpi, bbox_inches='tight')
    print(f"Plot saved as total_fields.png (DPI: {dpi})")

    plt.show()

    print(f"\nDone! Peak memory: {get_memory_usage():.1f} MB")

if __name__ == "__main__":
    try:
        main()
    except MemoryError:
        print("\n" + "="*70)
        print("MEMORY ERROR!")
        print("="*70)
        print("The grid size is too large for available memory.")
        print("Solutions:")
        print("  1. Close other applications to free RAM")
        print("  2. Reduce disc_per_lambda in mom.c (e.g., 40 → 30)")
        print("  3. The script will auto-adjust grid size on next run")
        print("="*70)
    except KeyboardInterrupt:
        print("\nInterrupted by user")
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()



