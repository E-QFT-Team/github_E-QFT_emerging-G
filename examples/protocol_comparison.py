#!/usr/bin/env python3
"""
Example comparing Protocol A (surface minimization) vs Protocol B (mass-defect).

This demonstrates both approaches to extracting the emergent gravitational
coupling G_eff from quantum field theory calculations.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add src and analysis to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))

from projector_factory import generate_fourier_projectors
from metric import set_projectors, commutator_metric, clear_cache


def protocol_a_simple(lattice_size=6):
    """
    Simplified Protocol A implementation for demonstration.
    
    Real implementation uses NetworkX min-cut - see analysis/final_comprehensive_validation.py
    """
    print(f"\n🔷 PROTOCOL A: Surface Minimization ({lattice_size}³)")
    print("-" * 50)
    
    # Generate projectors
    lattice_shape = (lattice_size,) * 3
    projectors = generate_fourier_projectors(lattice_shape, dimension=32)
    set_projectors(projectors)
    
    total_sites = lattice_size ** 3
    
    # Simple cubic region (half the lattice)
    region_sites = []
    for i in range(lattice_size):
        for j in range(lattice_size):
            for k in range(lattice_size // 2):  # Half in z-direction
                site_idx = i * lattice_size**2 + j * lattice_size + k
                region_sites.append(site_idx)
    
    print(f"   Region A: {len(region_sites)} sites")
    print(f"   Region B: {total_sites - len(region_sites)} sites")
    
    # Compute boundary entropy (simplified)
    boundary_pairs = []
    surface_entropy = 0
    
    for site_a in region_sites:
        # Check connections to complement
        for site_b in range(total_sites):
            if site_b not in region_sites:
                # Check if they're neighbors (simplified boundary detection)
                # Convert back to coordinates
                coords_a = np.unravel_index(site_a, lattice_shape)
                coords_b = np.unravel_index(site_b, lattice_shape)
                
                # Chebyshev distance
                distance = max(abs(coords_a[i] - coords_b[i]) for i in range(3))
                
                if distance == 1:  # Nearest neighbors across boundary
                    boundary_pairs.append((site_a, site_b))
                    
                    # Add contribution to entropy (simplified)
                    d_squared = commutator_metric(site_a, site_b)
                    surface_entropy += np.log(1 + d_squared)  # Simplified entropy
    
    print(f"   Boundary pairs: {len(boundary_pairs)}")
    print(f"   Surface entropy: {surface_entropy:.3f}")
    
    # Estimate surface area
    surface_area = len(boundary_pairs)  # Simplified
    
    # G_eff from area law
    if surface_area > 0:
        G_eff_A = surface_entropy / surface_area
    else:
        G_eff_A = 0
    
    print(f"   ✅ G_eff(A) = {G_eff_A:.4f} lattice units")
    
    return {
        'method': 'Protocol A',
        'G_eff': G_eff_A,
        'surface_entropy': surface_entropy,
        'surface_area': surface_area,
        'boundary_pairs': len(boundary_pairs)
    }


def protocol_b_simple(lattice_size=6):
    """
    Simplified Protocol B implementation for demonstration.
    
    Real implementation uses multi-site averaging - see analysis/protocol_b_stabilization.py
    """
    print(f"\n🔶 PROTOCOL B: Mass-Defect Analysis ({lattice_size}³)")
    print("-" * 50)
    
    # Generate baseline projectors
    lattice_shape = (lattice_size,) * 3
    baseline_projectors = generate_fourier_projectors(lattice_shape, dimension=32)
    
    # Choose central site for perturbation
    center_idx = (lattice_size**3) // 2
    center_coords = np.unravel_index(center_idx, lattice_shape)
    print(f"   Central site: {center_idx} at {center_coords}")
    
    # Test different μ values
    mu_values = [0.05, 0.1, 0.2]
    kappa_values = []
    
    for mu in mu_values:
        print(f"\n   Testing μ = {mu}")
        
        # Create perturbed projectors
        perturbed_projectors = baseline_projectors.copy()
        perturbed_projectors[center_idx] = (1 + mu) * baseline_projectors[center_idx]
        
        # Compute h₀₀(r) profile
        set_projectors(baseline_projectors)
        baseline_distances = {}
        
        set_projectors(perturbed_projectors)
        perturbed_distances = {}
        
        # Collect distances at different radii
        for radius in [1, 2, 3]:
            baseline_list = []
            perturbed_list = []
            
            for site in range(lattice_size**3):
                site_coords = np.unravel_index(site, lattice_shape)
                site_radius = max(abs(site_coords[i] - center_coords[i]) for i in range(3))
                
                if site_radius == radius:
                    # Baseline
                    set_projectors(baseline_projectors)
                    d_base = commutator_metric(center_idx, site)
                    baseline_list.append(d_base)
                    
                    # Perturbed  
                    set_projectors(perturbed_projectors)
                    d_pert = commutator_metric(center_idx, site)
                    perturbed_list.append(d_pert)
            
            if baseline_list:
                baseline_distances[radius] = np.mean(baseline_list)
                perturbed_distances[radius] = np.mean(perturbed_list)
        
        # Compute h₀₀(r) = perturbed - baseline
        h00_profile = {}
        for radius in baseline_distances:
            h00_profile[radius] = perturbed_distances[radius] - baseline_distances[radius]
            print(f"      r={radius}: h₀₀ = {h00_profile[radius]:.6f}")
        
        # Fit κ from h₀₀(r) = -2κ/r
        if len(h00_profile) >= 2:
            radii = np.array(list(h00_profile.keys()))
            h00_vals = np.array(list(h00_profile.values()))
            
            # Simple linear fit: h₀₀ * r = -2κ
            kappa_estimates = -h00_vals * radii / 2
            kappa = np.mean(kappa_estimates)
            kappa_values.append(kappa)
            
            print(f"      ✅ κ = {kappa:.4f}")
        else:
            kappa_values.append(np.nan)
    
    # Check linearity κ ∝ μ
    valid_mu = []
    valid_kappa = []
    for mu, kappa in zip(mu_values, kappa_values):
        if not np.isnan(kappa):
            valid_mu.append(mu)
            valid_kappa.append(abs(kappa))
    
    if len(valid_mu) >= 2:
        # Linear fit
        coeffs = np.polyfit(valid_mu, valid_kappa, 1)
        slope = coeffs[0]
        intercept = coeffs[1]
        
        print(f"\n   Linear fit: |κ| = {slope:.3f} * μ + {intercept:.3f}")
        
        # Average κ as G_eff estimate
        G_eff_B = np.mean(valid_kappa)
    else:
        G_eff_B = np.nan
        slope = np.nan
    
    print(f"   ✅ G_eff(B) = {G_eff_B:.4f} lattice units")
    
    return {
        'method': 'Protocol B',
        'G_eff': G_eff_B,
        'mu_values': mu_values,
        'kappa_values': kappa_values,
        'slope': slope
    }


def compare_protocols():
    """Compare Protocol A vs Protocol B results."""
    
    print("🚀 PROTOCOL COMPARISON EXAMPLE")
    print("=" * 60)
    print("Comparing surface minimization vs mass-defect approaches")
    
    # Run both protocols
    result_A = protocol_a_simple(lattice_size=5)  # Smaller for demo
    result_B = protocol_b_simple(lattice_size=5)
    
    # Comparison
    print(f"\n📊 COMPARISON RESULTS")
    print("=" * 30)
    print(f"Protocol A G_eff: {result_A['G_eff']:.4f}")
    print(f"Protocol B G_eff: {result_B['G_eff']:.4f}")
    
    if not np.isnan(result_A['G_eff']) and not np.isnan(result_B['G_eff']):
        relative_diff = abs(result_A['G_eff'] - result_B['G_eff']) / result_A['G_eff'] * 100
        print(f"Relative difference: {relative_diff:.1f}%")
        
        if relative_diff < 20:
            print("✅ Protocols agree reasonably well")
        else:
            print("⚠️  Large discrepancy - may need refinement")
    
    # Create comparison plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: G_eff comparison
    methods = ['Protocol A\n(Surface)', 'Protocol B\n(Mass-defect)']
    G_values = [result_A['G_eff'], result_B['G_eff']]
    colors = ['blue', 'orange']
    
    bars = ax1.bar(methods, G_values, color=colors, alpha=0.7)
    ax1.set_ylabel('G_eff (lattice units)')
    ax1.set_title('Protocol Comparison')
    ax1.grid(True, alpha=0.3)
    
    # Add value labels
    for bar, value in zip(bars, G_values):
        if not np.isnan(value):
            ax1.text(bar.get_x() + bar.get_width()/2, value + 0.001, 
                    f'{value:.3f}', ha='center', va='bottom')
    
    # Plot 2: Protocol B linearity (if available)
    if not np.isnan(result_B['slope']):
        mu_vals = result_B['mu_values']
        kappa_vals = [abs(k) if not np.isnan(k) else 0 for k in result_B['kappa_values']]
        
        ax2.plot(mu_vals, kappa_vals, 'o-', color='orange', markersize=8, linewidth=2)
        ax2.set_xlabel('Mass-defect parameter μ')
        ax2.set_ylabel('|κ| parameter')
        ax2.set_title('Protocol B Linearity')
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'Protocol B\ndata unavailable', 
                ha='center', va='center', transform=ax2.transAxes)
    
    plt.tight_layout()
    plt.savefig('protocol_comparison_example.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"\n📈 Comparison plot saved as 'protocol_comparison_example.png'")
    print(f"\n📚 For full validation, run:")
    print(f"   python analysis/final_comprehensive_validation.py")
    
    return result_A, result_B


if __name__ == "__main__":
    results = compare_protocols()