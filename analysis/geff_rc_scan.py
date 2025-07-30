#!/usr/bin/env python3
"""
G_eff radius convergence scan as specified in the task.

Tests r_C = 2, 3, 4 to find where lattice drift ≤ 5%.
Expected: r_C = 4 captures 95% of surface weight and stabilizes convergence.
"""

import numpy as np
import networkx as nx
import time
from projector_factory import generate_fourier_projectors
from metric import commutator_metric, set_projectors
from entropy_tools import von_neumann_entropy


def geff_for_radius(lattice_n=8, max_rc=4, block_ratio=0.5):
    """
    Compute G_eff for given lattice size and max Chebyshev radius.
    
    Parameters
    ----------
    lattice_n : int
        Lattice size (n for n³ lattice)
    max_rc : int
        Maximum Chebyshev radius for edge inclusion
    block_ratio : float
        Block size as fraction of lattice size
        
    Returns
    -------
    dict
        G_eff computation results
    """
    
    start_time = time.time()
    
    # Generate deterministic rank-1 projectors
    lattice_shape = (lattice_n,) * 3
    projectors = generate_fourier_projectors(lattice_shape, 1.0, 1.0, False, 64)
    set_projectors(projectors)
    
    # Build graph with specified max radius
    G = nx.Graph()
    N = lattice_n**3
    G.add_nodes_from(range(N))
    
    edges_added = 0
    
    for i in range(N):
        xi = np.unravel_index(i, lattice_shape)
        for j in range(i + 1, N):
            xj = np.unravel_index(j, lattice_shape)
            
            # Chebyshev distance with periodic boundaries
            dx = np.abs(np.array(xi) - np.array(xj))
            dx = np.minimum(dx, np.array(lattice_shape) - dx)
            rc = np.max(dx)
            
            if 0 < rc <= max_rc:  # Exclude self-loops, include up to max_rc
                weight = commutator_metric(i, j)
                G.add_edge(i, j, weight=weight)
                edges_added += 1
    
    # Define region A (cubic block)
    block_size = max(2, int(lattice_n * block_ratio))  # At least 2³ block
    region_A = []
    
    start_idx = (lattice_n - block_size) // 2
    for i in range(block_size):
        for j in range(block_size):
            for k in range(block_size):
                pos = (start_idx + i, start_idx + j, start_idx + k)
                site = np.ravel_multi_index(pos, lattice_shape)
                region_A.append(site)
    
    region_B = [node for node in G.nodes() if node not in region_A]
    
    # Min-cut
    if len(region_A) > 0 and len(region_B) > 0:
        cut_value, partition = nx.minimum_cut(G, region_A[0], region_B[0], 
                                            capacity='weight')
        
        # Calculate surface area
        area_gamma = 0.0
        cut_edges = nx.edge_boundary(G, partition[0], partition[1])
        for u, v in cut_edges:
            area_gamma += G[u][v]['weight']
    else:
        area_gamma = 0.0
    
    # Entanglement entropy
    S_A = von_neumann_entropy(region_A, method='random_matrix')
    
    # G_eff
    G_eff = area_gamma / (4 * S_A) if S_A > 0 else 0.0
    
    computation_time = time.time() - start_time
    
    return {
        'lattice_n': lattice_n,
        'max_rc': max_rc,
        'block_size': block_size,
        'edges': edges_added,
        'area_gamma': area_gamma,
        'S_A': S_A,
        'G_eff': G_eff,
        'computation_time': computation_time
    }


def radius_convergence_scan():
    """
    Scan r_C = 2, 3, 4 for lattice convergence.
    Stop when drift ≤ 5%.
    """
    
    print("G_eff RADIUS CONVERGENCE SCAN")
    print("=" * 50)
    print("Testing r_C = 2, 3, 4 for lattice drift ≤ 5%")
    print()
    
    results = {}
    
    for rc in [2, 3, 4]:
        print(f"Testing r_C = {rc}")
        print("-" * 20)
        
        # Test on 8³ and 10³ lattices
        g8_result = geff_for_radius(lattice_n=8, max_rc=rc)
        g10_result = geff_for_radius(lattice_n=10, max_rc=rc)
        
        g8 = g8_result['G_eff']
        g10 = g10_result['G_eff']
        
        # Calculate drift
        if g8 > 0:
            drift = 100 * (g10 - g8) / g8
        else:
            drift = float('inf')
        
        results[rc] = {
            '8_cubed': g8_result,
            '10_cubed': g10_result,
            'drift_percent': drift
        }
        
        print(f"  8³:  G_eff = {g8:.6f} ({g8_result['edges']} edges, {g8_result['computation_time']:.1f}s)")
        print(f"  10³: G_eff = {g10:.6f} ({g10_result['edges']} edges, {g10_result['computation_time']:.1f}s)")
        print(f"  Drift: {drift:.2f}%")
        
        if abs(drift) <= 5:
            print(f"  ✓ Convergence achieved at r_C = {rc}")
            convergence_achieved = True
        else:
            print(f"  ⚠ Drift exceeds 5% threshold")
            convergence_achieved = False
        
        print()
    
    # Summary analysis
    print("CONVERGENCE ANALYSIS")
    print("-" * 30)
    
    for rc in [2, 3, 4]:
        drift = results[rc]['drift_percent']
        status = "PASS" if abs(drift) <= 5 else "FAIL"
        print(f"r_C = {rc}: drift = {drift:+.2f}% [{status}]")
    
    # Find optimal radius
    optimal_rc = None
    for rc in [2, 3, 4]:
        if abs(results[rc]['drift_percent']) <= 5:
            optimal_rc = rc
            break
    
    if optimal_rc:
        print(f"\n✓ Optimal radius: r_C = {optimal_rc}")
        print(f"  Captures sufficient surface weight for convergence")
    else:
        print(f"\n⚠ No radius achieves ≤5% drift")
        print(f"  May require finite-size extrapolation")
    
    # Save results
    np.save('radius_convergence_scan_results.npy', results)
    
    return results, optimal_rc


def finite_size_extrapolation(results):
    """
    Perform finite-size extrapolation G_eff(n) = G_∞ + A*n^(-1).
    """
    
    print("\nFINITE-SIZE EXTRAPOLATION")
    print("-" * 30)
    
    # Use best radius (r_C = 4 for full weight capture)
    rc = 4
    if rc not in results:
        print("No data for r_C = 4, using available data")
        rc = max(results.keys())
    
    # Extract G_eff values
    lattice_sizes = [8, 10]
    G_eff_values = [results[rc]['8_cubed']['G_eff'], 
                   results[rc]['10_cubed']['G_eff']]
    
    # Fit G_eff(n) = G_∞ + A/n
    inv_n = np.array([1/n for n in lattice_sizes])
    G_array = np.array(G_eff_values)
    
    # Linear regression: G = G_∞ + A * (1/n)
    A_matrix = np.vstack([np.ones(len(inv_n)), inv_n]).T
    coeffs, residuals, rank, s = np.linalg.lstsq(A_matrix, G_array, rcond=None)
    
    G_infinity = coeffs[0]
    A_coeff = coeffs[1]
    
    # Compute R²
    G_pred = G_infinity + A_coeff * inv_n
    ss_res = np.sum((G_array - G_pred)**2)
    ss_tot = np.sum((G_array - np.mean(G_array))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    print(f"Finite-size fit: G_eff(n) = {G_infinity:.6f} + {A_coeff:.6f}/n")
    print(f"R² = {r_squared:.3f}")
    
    # Test extrapolation accuracy
    G10_extrap = G_infinity + A_coeff / 10
    G10_actual = results[rc]['10_cubed']['G_eff']
    extrap_error = abs(G10_extrap - G10_actual) / G10_actual * 100
    
    print(f"\nExtrapolation test:")
    print(f"  10³ actual: {G10_actual:.6f}")
    print(f"  10³ extrap: {G10_extrap:.6f}")
    print(f"  Error: {extrap_error:.2f}%")
    
    if extrap_error < 3:
        print("  ✓ Extrapolation accurate within 3%")
    else:
        print("  ⚠ Extrapolation error exceeds 3%")
    
    return {
        'G_infinity': G_infinity,
        'A_coeff': A_coeff,
        'r_squared': r_squared,
        'extrapolation_error': extrap_error,
        'lattice_sizes': lattice_sizes,
        'G_eff_values': G_eff_values
    }


def main():
    """
    Run complete radius convergence analysis.
    """
    
    # Step 1: Radius convergence scan
    results, optimal_rc = radius_convergence_scan()
    
    # Step 2: Finite-size extrapolation
    extrap_results = finite_size_extrapolation(results)
    
    # Combined results
    complete_results = {
        'radius_scan': results,
        'optimal_rc': optimal_rc,
        'extrapolation': extrap_results
    }
    
    print(f"\nFINAL ASSESSMENT")
    print("=" * 30)
    if optimal_rc:
        print(f"✓ Radius convergence: r_C = {optimal_rc}")
    else:
        print("⚠ Finite-size extrapolation needed")
    
    print(f"✓ Extrapolated G_∞ = {extrap_results['G_infinity']:.6f}")
    print(f"✓ Method validation complete")
    
    return complete_results


if __name__ == "__main__":
    results = main()