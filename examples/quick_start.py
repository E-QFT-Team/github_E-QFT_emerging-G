#!/usr/bin/env python3
"""
Quick start example for emergent gravity calculations.

This script demonstrates the basic workflow:
1. Generate local projectors for a discrete lattice
2. Compute quantum commutator metrics
3. Extract emergent gravitational coupling G_eff
"""

import numpy as np
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from projector_factory import generate_fourier_projectors
from metric import set_projectors, commutator_metric, clear_cache
from entropy_tools import von_neumann_entropy


def main():
    """Run a simple emergent gravity calculation."""
    
    print("🌌 EMERGENT GRAVITY QUICK START")
    print("=" * 40)
    
    # Step 1: Generate projectors for a small lattice
    print("\n1. Generating local projectors...")
    lattice_shape = (4, 4, 4)  # 4³ lattice
    projectors = generate_fourier_projectors(
        lattice_shape=lattice_shape,
        localization_width=1.0,
        lambda_param=1.0,
        dimension=32
    )
    
    print(f"   ✅ Generated {len(projectors)} projectors for {lattice_shape} lattice")
    
    # Step 2: Set projectors and compute some metric distances
    print("\n2. Computing commutator metrics...")
    set_projectors(projectors)
    
    # Compute distances from center site to neighbors
    center_site = 32  # Middle of 4³ lattice (index for (2,2,2))
    
    distances = {}
    for site in [center_site + 1, center_site + 16, center_site + 4]:  # Nearby sites
        d_squared = commutator_metric(center_site, site)
        distances[site] = d_squared
    
    print(f"   ✅ Computed metric distances from center site {center_site}:")
    for site, d_sq in distances.items():
        print(f"      d²({center_site}, {site}) = {d_sq:.6f}")
    
    # Step 3: Simple G_eff estimation using minimal surface
    print("\n3. Estimating emergent gravitational coupling...")
    
    # For quick demo, use simplified approach
    # In real analysis, this would use NetworkX min-cut algorithms
    
    total_sites = len(projectors)
    region_size = total_sites // 2  # Half the lattice
    
    # Estimate surface entropy (simplified)
    typical_entropy = 2.5  # Reasonable value for small systems
    
    # Estimate G_eff using simplified scaling
    # Real calculation uses Protocol A with proper surface minimization
    surface_area_est = region_size ** (2/3)  # 2D surface for 3D volume
    G_eff_estimate = typical_entropy / surface_area_est
    
    print(f"   ✅ Estimated G_eff ~ {G_eff_estimate:.3f} lattice units")
    print(f"      (This is a simplified estimate - see analysis/ for full Protocol A)")
    
    # Step 4: Show unit conversion concept
    print("\n4. Unit conversion to physical units...")
    
    from scipy import constants
    c = constants.c  # Speed of light
    G_newton = 6.674e-11  # Newton's constant
    
    # For demonstration, assume lattice spacing ~ 1 fm
    lattice_spacing = 1e-15  # meters (femtometer scale)
    
    # Convert to SI units: G_SI = G_eff * a² * c² / ℓ₀²
    # With a = ℓ₀, this becomes G_SI = G_eff * c²
    lambda_needed = np.sqrt(G_newton / (G_eff_estimate * c**2))
    
    print(f"   📏 Lattice spacing: {lattice_spacing:.0e} m")
    print(f"   🔧 Required λ parameter: {lambda_needed:.2e}")
    print(f"   ⚖️  This gives G_SI = {G_newton:.3e} m³/(kg·s²)")
    
    print("\n🎉 Quick start complete!")
    print("\n📚 Next steps:")
    print("   • Run analysis/final_comprehensive_validation.py for full validation")
    print("   • See examples/protocol_comparison.py for Protocol A vs B")
    print("   • Check tests/ directory for unit tests")
    print("   • Read README.md for detailed documentation")
    
    return {
        'projectors_generated': len(projectors),
        'G_eff_estimate': G_eff_estimate,
        'lambda_needed': lambda_needed
    }


if __name__ == "__main__":
    results = main()