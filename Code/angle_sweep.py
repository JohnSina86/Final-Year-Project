"""
L-Bracket MoM Angle Sweep Orchestrator - OPTIMIZED VERSION
Runs parametric study of optimal omega vs geometry angle
"""

import os
import sys
import subprocess
import time
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
from multiprocessing import Pool, cpu_count
import re

# Progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    print("Note: Install tqdm for better progress bars: pip install tqdm")

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIG = {
    'angle_start': 0,
    'angle_end': 63,
    'angle_step': 0.1,  # TEST MODE: 10° increments
    'executable': './mom_angle_sweep.exe',
    'num_cores': 8,
    'parallel_jobs': 1,  # Run 8 angles at once (1 core each)
    'results_dir': 'results',
    'output_csv': 'omega_summary.csv',
    'plot_file': 'omega_vs_angle.png',
    'geometry_grid_file': 'geometry_grid.png',
    'save_full_output': False,  # Set True to keep all current_distribution files
    'visualize_geometry': True,  # Generate geometry plots
}

# ============================================================================
# GEOMETRY VISUALIZATION
# ============================================================================

def plot_geometry_from_file(geom_file, angle, output_file):
    """Read geometry.txt and create visualization"""
    try:
        # Read geometry file
        with open(geom_file, 'r') as f:
            lines = f.readlines()
        
        # Parse coordinates (skip header lines starting with %)
        x_coords = []
        y_coords = []
        for line in lines:
            if line.strip() and not line.startswith('%'):
                try:
                    x, y = map(float, line.split())
                    x_coords.append(x)
                    y_coords.append(y)
                except:
                    continue
        
        if len(x_coords) < 2:
            return False
        
        # Create plot
        plt.figure(figsize=(8, 8))
        plt.plot(x_coords, y_coords, 'b-', linewidth=2, label='L-Bracket')
        plt.plot(x_coords, y_coords, 'ro', markersize=3, alpha=0.5)
        
        # Mark start and end
        plt.plot(x_coords[0], y_coords[0], 'go', markersize=10, label='Start', zorder=5)
        plt.plot(x_coords[-1], y_coords[-1], 'rs', markersize=10, label='End', zorder=5)
        
        plt.xlabel('X Position (m)', fontsize=12)
        plt.ylabel('Y Position (m)', fontsize=12)
        plt.title(f'L-Bracket Geometry - Angle = {angle:.1f}°', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.axis('equal')
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        return True
        
    except Exception as e:
        print(f"    Warning: Could not plot geometry for angle {angle}: {e}")
        return False

def create_geometry_grid(df, results_dir):
    """Create grid of geometry visualizations"""
    print("="*70)
    print("  CREATING GEOMETRY VISUALIZATIONS")
    print("="*70)
    
    converged = df[df['converged'] == True].copy()
    
    if len(converged) == 0:
        print("  No successful simulations to visualize")
        return
    
    # Select up to 10 evenly spaced angles for visualization
    n_plots = min(10, len(converged))
    indices = np.linspace(0, len(converged)-1, n_plots, dtype=int)
    selected = converged.iloc[indices]
    
    # Calculate grid dimensions
    n_cols = min(5, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_plots == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    print(f"  Creating geometry grid with {n_plots} angles...")
    
    for idx, (_, row) in enumerate(selected.iterrows()):
        angle = row['angle']
        angle_dir = Path(results_dir) / f"angle_{int(angle*10):04d}"
        geom_file = angle_dir / 'geometry.txt'
        
        ax = axes[idx]
        
        if geom_file.exists():
            # Read and plot geometry
            try:
                with open(geom_file, 'r') as f:
                    lines = f.readlines()
                
                x_coords = []
                y_coords = []
                for line in lines:
                    if line.strip() and not line.startswith('%'):
                        try:
                            x, y = map(float, line.split())
                            x_coords.append(x)
                            y_coords.append(y)
                        except:
                            continue
                
                if len(x_coords) >= 2:
                    ax.plot(x_coords, y_coords, 'b-', linewidth=1.5)
                    ax.plot(x_coords, y_coords, 'ro', markersize=2, alpha=0.5)
                    ax.plot(x_coords[0], y_coords[0], 'go', markersize=6, zorder=5)
                    ax.set_title(f'{angle:.1f}° (ω={row["omega"]:.3f})', fontsize=10, fontweight='bold')
                    ax.set_xlabel('X (m)', fontsize=8)
                    ax.set_ylabel('Y (m)', fontsize=8)
                    ax.grid(True, alpha=0.3)
                    ax.set_aspect('equal')
                else:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'{angle:.1f}° - Error', fontsize=10)
            
            except Exception as e:
                ax.text(0.5, 0.5, f'Error:\n{str(e)[:30]}', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{angle:.1f}° - Error', fontsize=10)
        else:
            ax.text(0.5, 0.5, 'No geometry file', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{angle:.1f}° - Missing', fontsize=10)
    
    # Hide unused subplots
    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('L-Bracket Geometry Evolution with Angle', fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(CONFIG['geometry_grid_file'], dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Geometry grid saved: {CONFIG['geometry_grid_file']}")
    print()

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def print_header():
    """Print program header"""
    print("\n" + "="*70)
    print("  L-BRACKET ANGLE SWEEP ORCHESTRATOR - OPTIMIZED")
    print("  MoM Solver Parametric Study - PARALLEL MODE")
    print("="*70)
    print(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if CONFIG['angle_step'] >= 5:
        print(f"  MODE: TEST RUN (Step = {CONFIG['angle_step']}°)")
    print()

def run_simulation(angle):
    """Run single simulation for given angle - OPTIMIZED"""
    # Create unique working directory
    angle_dir = Path(CONFIG['results_dir']) / f"angle_{int(angle*10):04d}"
    angle_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy Transmitter pos file to subdirectory
    if os.path.exists('Transmitter pos'):
        shutil.copy('Transmitter pos', angle_dir / 'Transmitter pos')
    
    # Write angle file in subdirectory
    angle_file = angle_dir / 'angle.txt'
    with open(angle_file, 'w') as f:
        f.write(f"{angle}\n")
    
    # Set environment for OpenMP
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(max(1, CONFIG['num_cores'] // CONFIG['parallel_jobs']))
    
    start_time = time.time()
    try:
        # Run in subdirectory to avoid conflicts
        result = subprocess.run(
            [os.path.abspath(CONFIG['executable'])],
            capture_output=True,
            text=True,
            timeout=1000,
            cwd=str(angle_dir),
            env=env
        )
        elapsed = time.time() - start_time
        
        if result.returncode != 0:
            # Save stderr for debugging
            with open(angle_dir / 'error.log', 'w') as f:
                f.write(result.stderr)
            return {
                'angle': angle,
                'omega': None,
                'iterations': None,
                'time': elapsed,
                'error': None,
                'converged': False,
                'status': f'ERROR: Return code {result.returncode}'
            }
        
        # Fast parsing with regex
        omega, iterations, error = parse_output_fast(result.stdout)
        
        # Create individual geometry plot if requested
        if CONFIG['visualize_geometry']:
            geom_file = angle_dir / 'geometry.txt'
            if geom_file.exists():
                plot_file = angle_dir / f'geometry_{int(angle*10):04d}.png'
                plot_geometry_from_file(geom_file, angle, plot_file)
        
        # Optionally delete large text files to save disk space
        if not CONFIG['save_full_output']:
            for f in ['current_distribution.txt']:
                fpath = angle_dir / f
                if fpath.exists():
                    fpath.unlink()
        
        return {
            'angle': angle,
            'omega': omega,
            'iterations': iterations,
            'time': elapsed,
            'error': error,
            'converged': True,
            'status': 'SUCCESS'
        }
        
    except subprocess.TimeoutExpired:
        return {
            'angle': angle,
            'omega': None,
            'iterations': None,
            'time': 300,
            'error': None,
            'converged': False,
            'status': 'ERROR: Timeout'
        }
    except Exception as e:
        return {
            'angle': angle,
            'omega': None,
            'iterations': None,
            'time': time.time() - start_time,
            'error': None,
            'converged': False,
            'status': f'ERROR: {str(e)}'
        }

def parse_output_fast(stdout):
    """Fast regex-based parsing"""
    omega_match = re.search(r'OPTIMAL OMEGA:\s+([\d.]+)', stdout)
    omega = float(omega_match.group(1)) if omega_match else None
    
    iter_match = re.search(r'Converged in (\d+) iterations', stdout)
    iterations = int(iter_match.group(1)) if iter_match else None
    
    error_match = re.search(r'error=([\d.e+-]+)', stdout)
    error = float(error_match.group(1)) if error_match else None
    
    return omega, iterations, error

def run_all_angles_parallel():
    """Run simulations in parallel - OPTIMIZED"""
    
    # Generate angle list
    angles = []
    angle = CONFIG['angle_start']
    while angle <= CONFIG['angle_end']:
        angles.append(round(angle, 2))
        angle += CONFIG['angle_step']
    
    print("="*70)
    print("  STARTING PARALLEL SIMULATIONS")
    print("="*70)
    print(f"  Angles: {CONFIG['angle_start']}° to {CONFIG['angle_end']}° (step={CONFIG['angle_step']}°)")
    print(f"  Total simulations: {len(angles)}")
    print(f"  Parallel jobs: {CONFIG['parallel_jobs']} ({max(1, CONFIG['num_cores']//CONFIG['parallel_jobs'])} cores each)")
    print(f"  Executable: {CONFIG['executable']}")
    print()
    
    # Check executable exists
    if not os.path.exists(CONFIG['executable']):
        print(f"  ERROR: Executable '{CONFIG['executable']}' not found!")
        print("  Please compile your C code first:")
        print("    gcc -std=c11 -O3 -march=native -fopenmp -ffast-math \\")
        print("        -funroll-loops mom_angle_sweep.c -lm -o mom_angle_sweep.exe")
        return None
    
    # Check config file exists
    if not os.path.exists('Transmitter pos'):
        print("  ERROR: 'Transmitter pos' file not found!")
        return None
    
    # Create results directory
    Path(CONFIG['results_dir']).mkdir(exist_ok=True)
    
    # Run in parallel
    start_time = time.time()
    
    print("  Running simulations...\n")
    
    with Pool(processes=CONFIG['parallel_jobs']) as pool:
        results = []
        converged_count = 0
        failed_count = 0
        
        # Use tqdm if available, otherwise fallback
        if HAS_TQDM:
            pbar = tqdm(total=len(angles), 
                       desc="  Progress",
                       unit="sim",
                       ncols=90,
                       mininterval=0.01,
                       maxinterval=0.1,
                       miniters=1,
                       smoothing=0.1,
                       bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}')
            
            for result in pool.imap_unordered(run_simulation, angles):
                results.append(result)
                if result['converged']:
                    converged_count += 1
                    omega_val = result['omega']
                else:
                    failed_count += 1
                    omega_val = None
                
                pbar.set_postfix({
                    '✓': converged_count, 
                    '✗': failed_count,
                    'ω': f"{omega_val:.3f}" if omega_val else "N/A"
                }, refresh=True)
                pbar.update(1)
            
            pbar.close()
        else:
            # Fallback without tqdm
            completed = 0
            last_update = 0
            update_interval = max(1, len(angles) // 100)
            
            for result in pool.imap_unordered(run_simulation, angles):
                results.append(result)
                completed += 1
                
                if result['converged']:
                    converged_count += 1
                else:
                    failed_count += 1
                
                if completed - last_update >= update_interval or len(angles) < 100:
                    pct = 100 * completed / len(angles)
                    elapsed = time.time() - start_time
                    rate = completed / elapsed if elapsed > 0 else 0
                    eta = (len(angles) - completed) / rate if rate > 0 else 0
                    
                    bar_length = 40
                    filled = int(bar_length * completed / len(angles))
                    bar = '█' * filled + '░' * (bar_length - filled)
                    
                    print(f"\r  [{bar}] {pct:5.1f}% | {completed}/{len(angles)} | "
                          f"✓{converged_count} ✗{failed_count} | "
                          f"{rate:.1f}/s | ETA:{eta/60:.1f}m   ", end='', flush=True)
                    
                    last_update = completed
            
            print()
    
    total_time = time.time() - start_time
    
    print()
    print("="*70)
    print(f"  SWEEP COMPLETE")
    print("="*70)
    print(f"  Total execution time: {int(total_time//60)}m {int(total_time%60)}s")
    print(f"  Average rate: {len(angles)/total_time:.2f} simulations/second")
    print(f"  Speedup vs sequential: ~{CONFIG['parallel_jobs']}x")
    print(f"  Success rate: {converged_count}/{len(angles)} ({100*converged_count/len(angles):.1f}%)")
    print()
    
    # Sort by angle before returning
    df = pd.DataFrame(results)
    df = df.sort_values('angle').reset_index(drop=True)
    return df

def save_results(df):
    """Save results to CSV"""
    print("="*70)
    print("  RESULTS SAVED")
    print("="*70)
    print(f"  File: {CONFIG['output_csv']}")
    print(f"  Rows: {len(df)}")
    print()
    
    df.to_csv(CONFIG['output_csv'], index=False)
    
    # Show sample
    print("  Sample (first 10 rows):")
    print(df.head(10).to_string(index=False))
    print()

def print_statistics(df):
    """Print summary statistics"""
    print("="*70)
    print("  SUMMARY STATISTICS")
    print("="*70)
    
    converged = df[df['converged'] == True]
    failed = df[df['converged'] == False]
    
    print(f"  Total simulations: {len(df)}")
    print(f"  Converged: {len(converged)} ({100*len(converged)/len(df):.1f}%)")
    print(f"  Failed: {len(failed)}")
    print()
    
    if len(converged) > 0:
        print("  Omega Statistics:")
        print(f"    Min:     {converged['omega'].min():.6f}")
        print(f"    Max:     {converged['omega'].max():.6f}")
        print(f"    Mean:    {converged['omega'].mean():.6f}")
        print(f"    Std Dev: {converged['omega'].std():.6f}")
        print(f"    Range:   {converged['omega'].max() - converged['omega'].min():.6f}")
        print()
        
        print("  Iterations Statistics:")
        print(f"    Min:     {int(converged['iterations'].min())}")
        print(f"    Max:     {int(converged['iterations'].max())}")
        print(f"    Mean:    {int(converged['iterations'].mean())}")
        print()
        
        print("  Time Statistics:")
        print(f"    Per angle (avg): {converged['time'].mean():.2f} seconds")
        print(f"    Total compute:   {converged['time'].sum():.1f} seconds")
    
    if len(failed) > 0:
        print()
        print("  Failed Angles:")
        for _, row in failed.head(20).iterrows():
            print(f"    Angle {row['angle']:5.1f}°: {row['status']}")
        if len(failed) > 20:
            print(f"    ... and {len(failed)-20} more")
    
    print()

def create_plots(df):
    """Create visualization plots"""
    print("="*70)
    print("  CREATING VISUALIZATIONS")
    print("="*70)
    
    converged = df[df['converged'] == True]
    
    if len(converged) < 2:
        print("  ERROR: Not enough converged data to plot")
        return
    
    # Omega vs Angle plot
    plt.figure(figsize=(12, 6))
    plt.plot(converged['angle'], converged['omega'], 'o-', linewidth=1.5, markersize=6)
    plt.xlabel('Angle (degrees)', fontsize=12)
    plt.ylabel('Optimal Omega', fontsize=12)
    plt.title('Optimal Omega vs Geometry Angle', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(CONFIG['plot_file'], dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Omega plot saved: {CONFIG['plot_file']}")
    
    # Create geometry grid
    if CONFIG['visualize_geometry']:
        create_geometry_grid(df, CONFIG['results_dir'])

def print_footer():
    """Print completion message"""
    print("="*70)
    print("  COMPLETE")
    print("="*70)
    print(f"  Results directory: {CONFIG['results_dir']}/")
    print(f"  Summary CSV: {CONFIG['output_csv']}")
    print(f"  Omega plot: {CONFIG['plot_file']}")
    if CONFIG['visualize_geometry']:
        print(f"  Geometry grid: {CONFIG['geometry_grid_file']}")
    print("="*70)
    
    if CONFIG['angle_step'] >= 5:
        print("\n  ⚠️  TEST MODE COMPLETE")
        print("  To run high-resolution sweep, change CONFIG['angle_step'] to 0.1")
    print()

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main execution"""
    try:
        print_header()
        
        # Run sweep
        df = run_all_angles_parallel()
        
        if df is None:
            return 1
        
        # Save and analyze
        save_results(df)
        print_statistics(df)
        create_plots(df)
        print_footer()
        
        return 0
        
    except KeyboardInterrupt:
        print("\n\n*** INTERRUPTED BY USER ***\n")
        return 130
    except Exception as e:
        print(f"\n\n  FATAL ERROR: {str(e)}\n")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    exit_code = main()
    sys.exit(exit_code)
