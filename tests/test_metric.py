#!/usr/bin/env python3
"""
Unit tests for the quantum commutator metric module.

Tests the computation of d²(i,j) = Tr([Πi, Πj]†[Πi, Πj]) and related
metric validation functions.
"""

import numpy as np
import pytest
from src.metric import (
    commutator_metric,
    set_projectors,
    clear_cache,
    validate_projector
)
from src.projector_factory import generate_1d_test_projectors


class TestMetric:
    """Test suite for quantum commutator metric computation."""
    
    def test_metric_symmetry(self):
        """Test that d²(i,j) = d²(j,i)."""
        projectors = generate_1d_test_projectors(4, dimension=16)
        set_projectors(projectors)
        
        for i in range(4):
            for j in range(4):
                d_ij = commutator_metric(i, j)
                d_ji = commutator_metric(j, i)
                assert abs(d_ij - d_ji) < 1e-12, f"Asymmetry: d({i},{j})={d_ij}, d({j},{i})={d_ji}"
    
    def test_metric_positivity(self):
        """Test that d²(i,j) ≥ 0 for all i,j."""
        projectors = generate_1d_test_projectors(5, dimension=16)
        set_projectors(projectors)
        
        for i in range(5):
            for j in range(5):
                d_squared = commutator_metric(i, j)
                assert d_squared >= 0, f"Negative metric: d²({i},{j}) = {d_squared}"
    
    def test_metric_diagonal(self):
        """Test that d²(i,i) = 0."""
        projectors = generate_1d_test_projectors(5, dimension=16)
        set_projectors(projectors)
        
        for i in range(5):
            d_ii = commutator_metric(i, i)
            assert abs(d_ii) < 1e-12, f"Non-zero diagonal: d²({i},{i}) = {d_ii}"
    
    def test_flat_space_scaling(self):
        """Test exact flat space scaling d²(i,j) = (i-j)²."""
        projectors = generate_1d_test_projectors(6, dimension=32)
        set_projectors(projectors)
        
        max_error = 0
        for i in range(6):
            for j in range(6):
                expected = (i - j) ** 2
                actual = commutator_metric(i, j)
                error = abs(actual - expected)
                max_error = max(max_error, error)
                
                assert error < 1e-10, f"Scaling error at ({i},{j}): {actual} vs {expected}"
        
        print(f"Max flat space error: {max_error:.2e}")
    
    def test_caching_behavior(self):
        """Test that metric caching works correctly."""
        projectors = generate_1d_test_projectors(3, dimension=8)
        set_projectors(projectors)
        
        # First call - should compute
        d1 = commutator_metric(0, 1)
        
        # Second call - should use cache (same result)
        d2 = commutator_metric(0, 1)
        assert abs(d1 - d2) < 1e-15
        
        # Clear cache and recompute
        clear_cache()
        d3 = commutator_metric(0, 1)
        assert abs(d1 - d3) < 1e-15
    
    def test_projector_validation(self):
        """Test projector validation function."""
        # Valid projector (Hermitian)
        P_valid = np.array([[1, 0], [0, 0]], dtype=complex)
        assert validate_projector(P_valid), "Valid projector rejected"
        
        # Invalid projector (non-Hermitian)
        P_invalid = np.array([[1, 1], [0, 0]], dtype=complex)
        assert not validate_projector(P_invalid), "Invalid projector accepted"
        
        # Nearly Hermitian (within tolerance)
        P_nearly = np.array([[1, 1e-14], [1e-14, 0]], dtype=complex)
        assert validate_projector(P_nearly), "Nearly Hermitian projector rejected"
    
    def test_large_projector_handling(self):
        """Test handling of larger projector matrices."""
        # Generate larger projectors
        from src.projector_factory import generate_fourier_projectors
        
        projectors = generate_fourier_projectors((3, 3, 3), dimension=64)
        set_projectors(projectors)
        
        # Test a few metric computations
        d_01 = commutator_metric(0, 1)
        d_13 = commutator_metric(1, 13)
        
        assert d_01 >= 0
        assert d_13 >= 0
        assert not np.isnan(d_01)
        assert not np.isnan(d_13)
    
    def test_metric_triangle_inequality(self):
        """Test approximate triangle inequality for metric."""
        projectors = generate_1d_test_projectors(4, dimension=16)
        set_projectors(projectors)
        
        # For sites 0, 1, 2: d(0,2) ≤ d(0,1) + d(1,2) approximately
        d_02 = np.sqrt(commutator_metric(0, 2))
        d_01 = np.sqrt(commutator_metric(0, 1))
        d_12 = np.sqrt(commutator_metric(1, 2))
        
        # Should be exactly d(0,2) = 2, d(0,1) = 1, d(1,2) = 1
        # So 2 ≤ 1 + 1 is satisfied
        assert d_02 <= d_01 + d_12 + 1e-10, "Triangle inequality violated"


if __name__ == "__main__":
    # Run basic tests
    test_suite = TestMetric()
    
    print("Running metric computation tests...")
    
    test_suite.test_metric_symmetry()
    print("✅ Metric symmetry")
    
    test_suite.test_metric_positivity()
    print("✅ Metric positivity")
    
    test_suite.test_metric_diagonal()
    print("✅ Zero diagonal")
    
    test_suite.test_flat_space_scaling()
    print("✅ Flat space scaling")
    
    test_suite.test_caching_behavior()
    print("✅ Caching behavior")
    
    test_suite.test_projector_validation()
    print("✅ Projector validation")
    
    test_suite.test_large_projector_handling()
    print("✅ Large projector handling")
    
    test_suite.test_metric_triangle_inequality()
    print("✅ Triangle inequality")
    
    print("\n🎉 All metric tests passed!")