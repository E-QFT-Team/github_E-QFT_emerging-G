#!/usr/bin/env python3
"""
Metric normalization as specified in task step 3.

Set d̄²_nn = 1 by normalizing all link weights by the mean nearest-neighbor distance.
This fixes G_eff to O(1) and makes λ-tuning transparent.

Formula: d²_ij ← d²_ij / d̄²_nn where d̄²_nn = (1/N_nn) Σ_<i,j> d²_ij
"""

import numpy as np
import networkx as nx
from projector_factory import generate_fourier_projectors
from metric import commutator_metric, set_projectors
from entropy_tools import von_neumann_entropy
import copy


def compute_mean_nn_distance(lattice_size=8):
    """
    Compute mean nearest-neighbor distance d̄²_nn.
    
    Parameters
    ----------
    lattice_size : int
        Lattice size (n for n³ lattice)
        
    Returns
    -------
    float
        Mean nearest-neighbor distance squared
    """
    
    print(f"Computing mean nearest-neighbor distance for {lattice_size}³ lattice")
    
    # Generate projectors
    lattice_shape = (lattice_size,) * 3
    projectors = generate_fourier_projectors(lattice_shape, 1.0, 1.0, False, 64)
    set_projectors(projectors)
    
    # Collect all nearest-neighbor distances
    nn_distances = []
    N = lattice_size**3
    
    for i in range(N):
        xi = np.unravel_index(i, lattice_shape)
        
        for j in range(i + 1, N):
            xj = np.unravel_index(j, lattice_shape)
            
            # Check if nearest neighbors (Chebyshev distance = 1)
            dx = np.abs(np.array(xi) - np.array(xj))
            dx = np.minimum(dx, np.array(lattice_shape) - dx)  # Periodic boundaries
            rc = np.max(dx)
            
            if rc == 1:  # Nearest neighbors only
                d_squared = commutator_metric(i, j)
                nn_distances.append(d_squared)
    
    d_bar_nn_squared = np.mean(nn_distances)
    
    print(f"  Found {len(nn_distances)} nearest-neighbor pairs")
    print(f"  Mean NN distance: d̄²_nn = {d_bar_nn_squared:.6f}")
    
    return d_bar_nn_squared, nn_distances


class NormalizedMetric:
    """
    Wrapper for metric computation with normalization.
    
    This class applies the normalization d²_ij ← d²_ij / d̄²_nn
    to all distance calculations.
    """
    
    def __init__(self, lattice_size=8):
        """
        Initialize normalized metric.
        
        Parameters
        ----------
        lattice_size : int
            Lattice size for computing normalization factor
        """
        self.d_bar_nn_squared, _ = compute_mean_nn_distance(lattice_size)
        self.lattice_size = lattice_size
        
        print(f"Normalized metric initialized: d̄²_nn = {self.d_bar_nn_squared:.6f}")
    
    def normalized_commutator_metric(self, i, j):
        """
        Compute normalized commutator metric.
        
        Parameters
        ----------
        i, j : int
            Site indices
            
        Returns
        -------
        float
            Normalized distance squared: d²_ij / d̄²_nn
        """
        raw_distance = commutator_metric(i, j)
        return raw_distance / self.d_bar_nn_squared


def geff_with_normalized_metric(lattice_n=8, max_rc=4, block_ratio=0.5):
    """
    Compute G_eff with normalized metric.
    
    This should yield G_eff ~ O(1) due to normalization.
    
    Parameters
    ----------
    lattice_n : int
        Lattice size
    max_rc : int
        Maximum Chebyshev radius
    block_ratio : float
        Block size ratio
        
    Returns
    -------
    dict
        G_eff results with normalized metric
    """
    
    print(f"Computing G_eff with normalized metric for {lattice_n}³ lattice")
    
    # Initialize normalized metric
    norm_metric = NormalizedMetric(lattice_n)
    
    # Generate projectors
    lattice_shape = (lattice_n,) * 3
    projectors = generate_fourier_projectors(lattice_shape, 1.0, 1.0, False, 64)
    set_projectors(projectors)
    
    # Build graph with normalized weights
    G = nx.Graph()
    N = lattice_n**3
    G.add_nodes_from(range(N))
    
    edges_added = 0
    
    for i in range(N):
        xi = np.unravel_index(i, lattice_shape)
        for j in range(i + 1, N):
            xj = np.unravel_index(j, lattice_shape)
            
            # Chebyshev distance
            dx = np.abs(np.array(xi) - np.array(xj))
            dx = np.minimum(dx, np.array(lattice_shape) - dx)
            rc = np.max(dx)
            
            if 0 < rc <= max_rc:
                # Use normalized metric
                weight = norm_metric.normalized_commutator_metric(i, j)
                G.add_edge(i, j, weight=weight)
                edges_added += 1
    
    print(f"  Graph: {edges_added} edges with normalized weights")
    
    # Define region A
    block_size = max(2, int(lattice_n * block_ratio))
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
        
        # Calculate normalized surface area
        area_gamma = 0.0
        cut_edges = nx.edge_boundary(G, partition[0], partition[1])
        for u, v in cut_edges:
            area_gamma += G[u][v]['weight']
    else:
        area_gamma = 0.0
    
    # Entanglement entropy
    S_A = von_neumann_entropy(region_A, method='random_matrix')
    
    # G_eff (now normalized)
    G_eff_normalized = area_gamma / (4 * S_A) if S_A > 0 else 0.0
    
    print(f"  Normalized Area(γ_A): {area_gamma:.6f}")
    print(f"  S(A): {S_A:.6f}")
    print(f"  G_eff (normalized): {G_eff_normalized:.6f}")
    
    return {
        'lattice_n': lattice_n,
        'max_rc': max_rc,
        'block_size': block_size,
        'normalization_factor': norm_metric.d_bar_nn_squared,
        'area_gamma_normalized': area_gamma,
        'S_A': S_A,
        'G_eff_normalized': G_eff_normalized,
        'edges': edges_added
    }


def protocol_b_with_normalized_metric(lattice_size=10, mu=0.2):
    """
    Run Protocol B with normalized metric.
    
    This should yield κ in the target range 1.8-2.2.
    
    Parameters
    ----------
    lattice_size : int
        Lattice size
    mu : float
        Mass-defect parameter
        
    Returns
    -------
    dict
        Protocol B results with normalized metric
    """
    
    print(f"\nProtocol B with normalized metric")
    print(f"Lattice: {lattice_size}³, μ = {mu}")
    
    # Initialize normalized metric
    norm_metric = NormalizedMetric(lattice_size)
    
    # Generate projectors
    lattice_shape = (lattice_size,) * 3
    baseline_projectors = generate_fourier_projectors(lattice_shape, 1.0, 1.0, False, 64)
    
    # Body-centered sites
    central_sites = []
    for i in [4, 5]:
        for j in [4, 5]:
            for k in [4, 5]:
                coord = (i, j, k)
                site_idx = np.ravel_multi_index(coord, lattice_shape)
                central_sites.append((site_idx, coord))
    
    # Perturbed projectors
    perturbed_projectors = copy.deepcopy(baseline_projectors)
    for site_idx, coord in central_sites:
        perturbed_projectors[site_idx] = (1 + mu) * baseline_projectors[site_idx]
    
    print(f"Applied μ = {mu} to {len(central_sites)} body-centered sites")
    
    # Compute h₀₀(r) with normalized metric
    set_projectors(baseline_projectors)
    baseline_profile = {}
    
    for radius in range(1, 7):
        distances_at_r = []
        
        for site_idx, coord in central_sites:
            for i in range(lattice_size**3):
                site_coord = np.unravel_index(i, lattice_shape)
                dr = max(abs(site_coord[0] - coord[0]), 
                        abs(site_coord[1] - coord[1]), 
                        abs(site_coord[2] - coord[2]))
                
                if dr == radius:
                    # Use normalized metric
                    d_squared = norm_metric.normalized_commutator_metric(site_idx, i)
                    distances_at_r.append(d_squared)
        
        if distances_at_r:
            baseline_profile[radius] = np.mean(distances_at_r)
    
    # Perturbed profile
    set_projectors(perturbed_projectors)
    perturbed_profile = {}
    
    for radius in range(1, 7):
        distances_at_r = []
        
        for site_idx, coord in central_sites:
            for i in range(lattice_size**3):
                site_coord = np.unravel_index(i, lattice_shape)
                dr = max(abs(site_coord[0] - coord[0]), 
                        abs(site_coord[1] - coord[1]), 
                        abs(site_coord[2] - coord[2]))
                
                if dr == radius:
                    # Use normalized metric
                    d_squared = norm_metric.normalized_commutator_metric(site_idx, i)
                    distances_at_r.append(d_squared)
        
        if distances_at_r:
            perturbed_profile[radius] = np.mean(distances_at_r)
    
    # Compute h₀₀(r)
    h00_profile = {}
    for radius in baseline_profile:
        if radius in perturbed_profile:
            h00 = perturbed_profile[radius] - baseline_profile[radius]
            h00_profile[radius] = h00
            print(f"  r = {radius}: h₀₀ = {h00:.6f}")
    
    # Fit κ in range 3-6
    r_fit = []
    h00_fit = []
    
    for radius in range(3, 7):
        if radius in h00_profile:
            r_fit.append(radius)
            h00_fit.append(h00_profile[radius])
    
    r_fit = np.array(r_fit)
    h00_fit = np.array(h00_fit)
    
    if len(r_fit) >= 3:
        # Fit h₀₀(r) = -2κ/r
        from scipy import optimize
        
        def fit_function(r, kappa):
            return -2 * kappa / r
        
        try:
            popt, pcov = optimize.curve_fit(fit_function, r_fit, h00_fit)
            kappa_normalized = popt[0]
            kappa_err = np.sqrt(pcov[0, 0]) if pcov.size > 0 else 0
            
            # G_eff(B) with normalized metric
            G_eff_B_normalized = kappa_normalized / (2 * mu)
            
            print(f"  Fitted κ (normalized): {kappa_normalized:.6f} ± {kappa_err:.6f}")
            print(f"  G_eff(B) (normalized): {G_eff_B_normalized:.6f}")
            
            fit_success = True
            
        except Exception as e:
            print(f"  Fit failed: {e}")
            kappa_normalized = np.nan
            G_eff_B_normalized = np.nan
            fit_success = False
    else:
        print(f"  Insufficient data for fitting")
        kappa_normalized = np.nan
        G_eff_B_normalized = np.nan
        fit_success = False
    
    return {
        'lattice_size': lattice_size,
        'mu': mu,
        'normalization_factor': norm_metric.d_bar_nn_squared,
        'kappa_normalized': kappa_normalized,
        'G_eff_B_normalized': G_eff_B_normalized,
        'fit_success': fit_success,
        'h00_profile': h00_profile
    }


def test_normalization_effects():
    """
    Test the effects of metric normalization on both protocols.
    """
    
    print("METRIC NORMALIZATION TEST")
    print("=" * 40)
    print("Testing both Protocol A and B with normalized metric")
    print("Expected: G_eff ~ O(1), κ in range 1.8-2.2")
    print()
    
    # Test Protocol A with normalization
    print("Protocol A with normalized metric:")
    protocol_a_8 = geff_with_normalized_metric(lattice_n=8, max_rc=4)
    protocol_a_10 = geff_with_normalized_metric(lattice_n=10, max_rc=4)
    
    # Calculate convergence
    G8 = protocol_a_8['G_eff_normalized']
    G10 = protocol_a_10['G_eff_normalized']
    drift = abs(G10 - G8) / G8 * 100 if G8 > 0 else float('inf')
    
    print(f"\nProtocol A convergence:")
    print(f"  8³: G_eff = {G8:.6f}")
    print(f"  10³: G_eff = {G10:.6f}")
    print(f"  Drift: {drift:.2f}%")
    
    # Test Protocol B with normalization
    print(f"\nProtocol B with normalized metric:")
    protocol_b_result = protocol_b_with_normalized_metric(lattice_size=10, mu=0.2)
    
    # Assessment
    print(f"\nNORMALIZATION ASSESSMENT")
    print("-" * 30)
    
    # Check if G_eff is O(1)
    if 0.1 <= G8 <= 10 and 0.1 <= G10 <= 10:
        print("✓ Protocol A: G_eff ~ O(1)")
        a_order_pass = True
    else:
        print("⚠ Protocol A: G_eff not O(1)")
        a_order_pass = False
    
    # Check if κ is in target range
    if protocol_b_result['fit_success']:
        kappa = abs(protocol_b_result['kappa_normalized'])
        if 1.8 <= kappa <= 2.2:
            print(f"✓ Protocol B: κ = {kappa:.3f} in range [1.8, 2.2]")
            b_range_pass = True
        else:
            print(f"⚠ Protocol B: κ = {kappa:.3f} outside range [1.8, 2.2]")
            b_range_pass = False
    else:
        print("⚠ Protocol B: Fitting failed")
        b_range_pass = False
    
    # Check protocol agreement
    if protocol_b_result['fit_success']:
        G_B = protocol_b_result['G_eff_B_normalized']
        agreement = abs(G_B - G8) / G8 * 100 if G8 > 0 else float('inf')
        print(f"Protocol A vs B: {agreement:.1f}% difference")
        
        if agreement <= 10:
            print("✓ Protocol agreement within 10%")
            agreement_pass = True
        else:
            print("⚠ Protocol agreement exceeds 10%")
            agreement_pass = False
    else:
        agreement_pass = False
    
    # Overall assessment
    overall_pass = a_order_pass and b_range_pass and agreement_pass
    
    print(f"\nOverall normalization: {'SUCCESS' if overall_pass else 'NEEDS REFINEMENT'}")
    
    # Save results
    normalization_results = {
        'protocol_a_8': protocol_a_8,
        'protocol_a_10': protocol_a_10,
        'protocol_b': protocol_b_result,
        'convergence_drift': drift,
        'a_order_pass': a_order_pass,
        'b_range_pass': b_range_pass,
        'agreement_pass': agreement_pass,
        'overall_pass': overall_pass
    }
    
    np.save('metric_normalization_results.npy', normalization_results)
    
    return normalization_results


def main():
    """
    Run complete metric normalization analysis.
    """
    
    results = test_normalization_effects()
    
    return results


if __name__ == "__main__":
    results = main()