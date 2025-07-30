#!/usr/bin/env python3
"""
Edge profile analysis for emergent gravity.

This module analyzes the decay of link weights w(r_C) with Chebyshev distance
to determine the optimal edge cutoff for G_eff calculations.

The 555% sensitivity between 1-nn and 2-nn suggests that second-neighbor 
links carry significant weight and cannot be ignored.
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from projector_factory import generate_fourier_projectors
from metric import commutator_metric, set_projectors


def link_profile(lattice_shape=(8, 8, 8), max_r=4):
    """
    Compute mean link weight vs Chebyshev distance.
    
    Parameters
    ----------
    lattice_shape : tuple
        Lattice dimensions (Lx, Ly, Lz)
    max_r : int
        Maximum Chebyshev distance to analyze
        
    Returns
    -------
    dict
        {r_C: mean_weight} for each Chebyshev distance
    """
    
    print(f"Computing link profile for {lattice_shape} lattice, max_r={max_r}")
    
    # Generate projectors
    proj = generate_fourier_projectors(lattice_shape, 1.0, 1.0, False, 64)
    set_projectors(proj)
    
    # Group weights by Chebyshev distance
    w = defaultdict(list)
    N = np.prod(lattice_shape)
    
    # Progress tracking
    pairs_computed = 0
    total_pairs = N * (N - 1) // 2
    
    for i in range(N):
        xi = np.unravel_index(i, lattice_shape)
        
        for j in range(i + 1, N):
            xj = np.unravel_index(j, lattice_shape)
            
            # Compute Chebyshev distance with periodic boundaries
            dx = np.abs(np.array(xi) - np.array(xj))
            dx = np.minimum(dx, np.array(lattice_shape) - dx)  # Periodic wrap
            rc = np.max(dx)
            
            if rc <= max_r:
                weight = commutator_metric(i, j)
                w[rc].append(weight)
            
            pairs_computed += 1
            
            # Progress update
            if pairs_computed % 10000 == 0:
                progress = 100 * pairs_computed / total_pairs
                print(f"  Progress: {progress:.1f}% ({pairs_computed}/{total_pairs})")
    
    # Compute mean weights
    profile = {r: np.mean(weights) for r, weights in w.items()}
    
    print(f"Link profile computed:")
    for r in sorted(profile.keys()):
        print(f"  r_C = {r}: <w> = {profile[r]:.6f} (n = {len(w[r])})")
    
    return profile, w


def analyze_weight_decay(profile):
    """
    Analyze the decay pattern w(r_C) ~ r_C^(-alpha).
    
    Parameters
    ----------
    profile : dict
        {r_C: mean_weight} profile data
        
    Returns
    -------
    dict
        Fit results including alpha exponent
    """
    
    print("\nAnalyzing weight decay pattern")
    
    # Extract data
    r_values = np.array(sorted(profile.keys()))
    w_values = np.array([profile[r] for r in r_values])
    
    # Remove r=0 if present (infinite weight)
    if r_values[0] == 0:
        r_values = r_values[1:]
        w_values = w_values[1:]
    
    if len(r_values) < 2:
        print("⚠ Insufficient data for decay analysis")
        return {'alpha': np.nan, 'r_squared': 0}
    
    # Log-log fit: log(w) = log(A) - alpha * log(r)
    log_r = np.log(r_values)
    log_w = np.log(w_values)
    
    # Linear regression
    A_matrix = np.vstack([log_r, np.ones(len(log_r))]).T
    coeffs, residuals, rank, s = np.linalg.lstsq(A_matrix, log_w, rcond=None)
    
    alpha = -coeffs[0]  # Negative because w ~ r^(-alpha)
    log_A = coeffs[1]
    A = np.exp(log_A)
    
    # Compute R²
    w_pred = A * r_values**(-alpha)
    ss_res = np.sum((w_values - w_pred)**2)
    ss_tot = np.sum((w_values - np.mean(w_values))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    print(f"  Decay fit: w(r) = {A:.4f} * r^(-{alpha:.2f})")
    print(f"  R² = {r_squared:.3f}")
    
    return {
        'alpha': alpha,
        'A': A,
        'r_squared': r_squared,
        'r_values': r_values,
        'w_values': w_values,
        'w_pred': w_pred
    }


def find_cutoff_radius(profile, target_fraction=0.95):
    """
    Find radius that captures target fraction of total weight.
    
    Parameters
    ----------
    profile : dict
        {r_C: mean_weight} profile data
    target_fraction : float
        Target fraction of total weight (default: 0.95)
        
    Returns
    -------
    dict
        Cutoff analysis results
    """
    
    print(f"\nFinding cutoff radius for {target_fraction*100}% weight capture")
    
    # Sort by radius
    r_values = np.array(sorted(profile.keys()))
    w_values = np.array([profile[r] for r in r_values])
    
    # Estimate number of pairs at each radius (approximate)
    # For periodic lattice, this is complex but we approximate
    lattice_size = 8  # Assume 8³ for now
    n_pairs = []
    
    for r in r_values:
        if r == 0:
            n_pairs.append(0)  # No pairs at distance 0
        else:
            # Rough estimate: surface area of cube at distance r
            if r == 1:
                n_pairs.append(6 * lattice_size**2)  # Face neighbors
            elif r == 2:
                n_pairs.append(12 * lattice_size**2)  # Edge neighbors (approx)
            else:
                n_pairs.append(8 * lattice_size**2)  # Corner neighbors (approx)
    
    n_pairs = np.array(n_pairs)
    
    # Total weight per radius shell
    total_weight_per_shell = w_values * n_pairs
    
    # Cumulative weight
    cumulative_weight = np.cumsum(total_weight_per_shell)
    total_weight = cumulative_weight[-1]
    cumulative_fraction = cumulative_weight / total_weight
    
    # Find cutoff radius
    cutoff_idx = np.where(cumulative_fraction >= target_fraction)[0]
    
    if len(cutoff_idx) == 0:
        cutoff_radius = r_values[-1]
        achieved_fraction = cumulative_fraction[-1]
    else:
        cutoff_radius = r_values[cutoff_idx[0]]
        achieved_fraction = cumulative_fraction[cutoff_idx[0]]
    
    print(f"  Cutoff radius: r_C = {cutoff_radius}")
    print(f"  Achieved fraction: {achieved_fraction:.3f}")
    
    return {
        'cutoff_radius': cutoff_radius,
        'achieved_fraction': achieved_fraction,
        'r_values': r_values,
        'cumulative_fraction': cumulative_fraction,
        'total_weight_per_shell': total_weight_per_shell
    }


def plot_edge_profile(profile, decay_fit, cutoff_analysis, save_path="edge_profile.png"):
    """
    Create comprehensive edge profile plots.
    
    Parameters
    ----------
    profile : dict
        Link weight profile
    decay_fit : dict
        Decay analysis results
    cutoff_analysis : dict
        Cutoff analysis results
    save_path : str
        Output plot file path
    """
    
    print(f"\nGenerating edge profile plots")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Weight vs radius (linear)
    r_vals = sorted(profile.keys())
    w_vals = [profile[r] for r in r_vals]
    
    ax1.plot(r_vals, w_vals, 'bo-', markersize=8, linewidth=2, label='Measured')
    ax1.set_xlabel('Chebyshev distance $r_C$')
    ax1.set_ylabel('Mean link weight $w(r_C)$')
    ax1.set_title('Link Weight vs Distance')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Weight vs radius (log-log)
    if 'r_values' in decay_fit and len(decay_fit['r_values']) > 0:
        ax2.loglog(decay_fit['r_values'], decay_fit['w_values'], 'bo', 
                  markersize=8, label='Measured')
        ax2.loglog(decay_fit['r_values'], decay_fit['w_pred'], 'r-', 
                  linewidth=2, label=f'Fit: $r^{{-{decay_fit["alpha"]:.2f}}}$')
        ax2.set_xlabel('Chebyshev distance $r_C$')
        ax2.set_ylabel('Mean link weight $w(r_C)$')
        ax2.set_title(f'Log-Log Decay (R² = {decay_fit["r_squared"]:.3f})')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
    
    # Plot 3: Cumulative weight fraction
    if 'r_values' in cutoff_analysis:
        r_cum = cutoff_analysis['r_values']
        frac_cum = cutoff_analysis['cumulative_fraction']
        cutoff_r = cutoff_analysis['cutoff_radius']
        
        ax3.plot(r_cum, frac_cum, 'go-', markersize=6, linewidth=2)
        ax3.axhline(y=0.95, color='red', linestyle='--', 
                   label='95% target')
        ax3.axvline(x=cutoff_r, color='red', linestyle='--', 
                   label=f'Cutoff: $r_C$ = {cutoff_r}')
        ax3.set_xlabel('Chebyshev distance $r_C$')
        ax3.set_ylabel('Cumulative weight fraction')
        ax3.set_title('Weight Accumulation')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        ax3.set_ylim(0, 1.05)
    
    # Plot 4: Weight per shell
    if 'total_weight_per_shell' in cutoff_analysis:
        r_shell = cutoff_analysis['r_values']
        weight_shell = cutoff_analysis['total_weight_per_shell']
        
        ax4.bar(r_shell, weight_shell, alpha=0.7, color='skyblue')
        ax4.set_xlabel('Chebyshev distance $r_C$')
        ax4.set_ylabel('Total weight in shell')
        ax4.set_title('Weight Distribution by Shell')
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Edge profile plots saved to {save_path}")


def main():
    """
    Run complete edge profile analysis.
    """
    
    print("EDGE PROFILE ANALYSIS")
    print("=" * 50)
    
    # Step 1: Compute link profile
    profile, raw_weights = link_profile(lattice_shape=(8, 8, 8), max_r=4)
    
    # Step 2: Analyze decay pattern
    decay_fit = analyze_weight_decay(profile)
    
    # Step 3: Find cutoff radius
    cutoff_analysis = find_cutoff_radius(profile, target_fraction=0.95)
    
    # Step 4: Generate plots
    plot_edge_profile(profile, decay_fit, cutoff_analysis)
    
    # Step 5: Save results
    results = {
        'profile': profile,
        'raw_weights': dict(raw_weights),  # Convert defaultdict
        'decay_fit': decay_fit,
        'cutoff_analysis': cutoff_analysis
    }
    
    np.save('edge_profile_results.npy', results)
    
    # Summary
    print(f"\nEDGE PROFILE SUMMARY")
    print("-" * 30)
    print(f"Decay exponent α = {decay_fit['alpha']:.2f}")
    print(f"95% weight cutoff: r_C = {cutoff_analysis['cutoff_radius']}")
    
    if decay_fit['alpha'] <= 3:
        print("⚠ α ≤ 3: Must include extended connectivity!")
    else:
        print("✓ α > 3: Nearest-neighbor approximation valid")
    
    return results


if __name__ == "__main__":
    results = main()