#!/usr/bin/env python3
"""
Unit tests for the projector factory module.

Tests the generation of local projectors, including Fourier-localization
and the deterministic rank-1 construction that preserves locality.
"""

import numpy as np
import pytest
from src.projector_factory import (
    generate_fourier_projectors,
    generate_1d_test_projectors,
    _create_projector_matrix
)


class TestProjectorFactory:
    """Test suite for projector generation functions."""
    
    def test_fourier_projectors_basic(self):
        """Test basic Fourier projector generation."""
        lattice_shape = (4, 4, 4)
        projectors = generate_fourier_projectors(lattice_shape, dimension=16)
        
        # Check correct number of projectors
        assert len(projectors) == 64  # 4³
        
        # Check each projector is properly shaped
        for site_idx, proj in projectors.items():
            assert proj.shape == (16, 16)
            assert np.allclose(proj, proj.conj().T)  # Hermitian
    
    def test_projector_idempotency(self):
        """Test that projectors are approximately idempotent."""
        lattice_shape = (3, 3, 3)
        projectors = generate_fourier_projectors(lattice_shape, dimension=8)
        
        for site_idx, proj in projectors.items():
            # P² should ≈ P for rank-1 projectors
            proj_squared = proj @ proj
            idempotency_error = np.linalg.norm(proj_squared - proj, 'fro')
            assert idempotency_error < 1e-10, f"Site {site_idx}: non-idempotent"
    
    def test_1d_flat_space_test(self):
        """Test the 1D flat space construction."""
        num_sites = 5
        projectors = generate_1d_test_projectors(num_sites, dimension=32)
        
        # Test exact distance scaling
        from src.metric import set_projectors, commutator_metric, clear_cache
        
        clear_cache()
        set_projectors(projectors)
        
        # Check d²(i,j) = (i-j)² exactly
        for i in range(num_sites):
            for j in range(num_sites):
                expected = (i - j) ** 2
                actual = commutator_metric(i, j)
                assert abs(actual - expected) < 1e-12, f"Sites {i},{j}: {actual} ≠ {expected}"
    
    def test_sparse_projectors(self):
        """Test sparse matrix generation."""
        lattice_shape = (3, 3, 3)
        projectors = generate_fourier_projectors(
            lattice_shape, dimension=8, use_sparse=True
        )
        
        import scipy.sparse as sp
        for proj in projectors.values():
            assert sp.issparse(proj)
            assert proj.shape == (8, 8)
    
    def test_projector_normalization(self):
        """Test that projector states are normalized."""
        # Test the internal function directly
        fourier_weight = np.random.complex128((10, 10))
        proj = _create_projector_matrix(fourier_weight, 1.0, 8, False)
        
        # For rank-1 projector P = |ψ⟩⟨ψ|, trace(P) = ⟨ψ|ψ⟩ = 1
        trace = np.trace(proj)
        assert abs(trace - 1.0) < 1e-12, f"Trace = {trace}, expected 1"
    
    def test_locality_preservation(self):
        """Test that deterministic projectors preserve locality."""
        lattice_shape = (4, 4, 4)
        projectors = generate_fourier_projectors(lattice_shape, dimension=16)
        
        from src.metric import set_projectors, commutator_metric, clear_cache
        
        clear_cache()
        set_projectors(projectors)
        
        # Test that nearby sites have stronger coupling than distant ones
        center_site = 32  # Middle of 4³ lattice (2,2,2)
        
        # Compare nearest neighbor vs far corner
        nearest = 33  # (2,2,3) - distance 1
        far = 0      # (0,0,0) - distance sqrt(12) ≈ 3.46
        
        d_near = commutator_metric(center_site, nearest)
        d_far = commutator_metric(center_site, far)
        
        # Locality means nearby sites should have smaller metric distance
        assert d_near < d_far, f"Locality violated: d_near={d_near}, d_far={d_far}"


if __name__ == "__main__":
    # Run basic tests
    test_suite = TestProjectorFactory()
    
    print("Running projector factory tests...")
    
    test_suite.test_fourier_projectors_basic()
    print("✅ Basic Fourier projector generation")
    
    test_suite.test_projector_idempotency()
    print("✅ Projector idempotency")
    
    test_suite.test_1d_flat_space_test()
    print("✅ 1D flat space test")
    
    test_suite.test_sparse_projectors()
    print("✅ Sparse matrix generation")
    
    test_suite.test_projector_normalization()
    print("✅ Projector normalization")
    
    test_suite.test_locality_preservation()
    print("✅ Locality preservation")
    
    print("\n🎉 All projector factory tests passed!")