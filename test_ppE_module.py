#!/usr/bin/env python3
"""
Pytest tests for ppE waveform module.

Validates that phase modifications are applied correctly
and match analytical formulas to 1e-6 rad precision.
"""

import pytest
import numpy as np
from ppE_waveform_module import ppEWaveformGenerator


class TestppEWaveformModule:
    """Test suite for ppE waveform module."""
    
    def setup_method(self):
        """Setup test fixtures."""
        self.generator = ppEWaveformGenerator()
        self.tolerance = 1e-6  # 1 microrad precision requirement
    
    def test_phase_ppE_corrections_analytical(self):
        """Test ppE phase corrections match analytical formulas."""
        
        # Test parameters
        f_test = np.array([50.0, 100.0, 200.0])  # Hz
        M_total = 50.0  # M_sun
        delta_phi_1 = 0.01
        delta_phi_1p5 = 0.02
        
        # Compute corrections using module
        phase_corr = self.generator.phase_ppE_corrections(
            f_test, M_total, delta_phi_1, delta_phi_1p5
        )
        
        # Analytical calculation
        M_total_sec = M_total * self.generator.M_sun * self.generator.G / self.generator.c**3
        
        expected_1PN = delta_phi_1 * (np.pi * M_total_sec * f_test)**(-1)
        expected_1p5PN = delta_phi_1p5 * (np.pi * M_total_sec * f_test)**(-2/3)
        expected_total = expected_1PN + expected_1p5PN
        
        # Check precision
        max_error = np.max(np.abs(phase_corr - expected_total))
        
        assert max_error < self.tolerance, f"Phase error {max_error:.2e} exceeds tolerance {self.tolerance:.2e}"
    
    def test_phase_ppE_zero_corrections(self):
        """Test that zero corrections give zero phase."""
        
        f_test = np.logspace(1, 3, 50)  # 10-1000 Hz
        M_total = 30.0
        
        phase_corr = self.generator.phase_ppE_corrections(
            f_test, M_total, delta_phi_1=0.0, delta_phi_1p5=0.0
        )
        
        assert np.allclose(phase_corr, 0.0), "Zero corrections should give zero phase"
    
    def test_chirp_mass_calculation(self):
        """Test chirp mass calculation."""
        
        # Known test case
        m1, m2 = 36.0, 29.0  # GW150914-like
        expected_chirp = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
        
        computed_chirp = self.generator.compute_chirp_mass(m1, m2)
        
        assert abs(computed_chirp - expected_chirp) < 1e-10, "Chirp mass calculation error"
    
    def test_total_mass_calculation(self):
        """Test total mass calculation."""
        
        m1, m2 = 30.0, 25.0
        expected_total = m1 + m2
        
        computed_total = self.generator.compute_total_mass(m1, m2)
        
        assert computed_total == expected_total, "Total mass calculation error"
    
    def test_waveform_generation_dimensions(self):
        """Test that waveform generation returns correct dimensions."""
        
        f_array = np.linspace(20, 500, 100)
        m1, m2 = 30.0, 25.0
        
        h_plus, h_cross = self.generator.generate_ppE_waveform(
            f_array, m1, m2, delta_phi_1=0.01, delta_phi_1p5=0.02
        )
        
        assert len(h_plus) == len(f_array), "h_plus length mismatch"
        assert len(h_cross) == len(f_array), "h_cross length mismatch"
        assert np.all(np.isfinite(h_plus[f_array > 0])), "h_plus contains non-finite values"
        assert np.all(np.isfinite(h_cross[f_array > 0])), "h_cross contains non-finite values"
    
    def test_ppE_parameter_scaling(self):
        """Test that ppE corrections scale linearly with parameters."""
        
        f_test = np.array([100.0])  # Single frequency
        M_total = 40.0
        
        # Test delta_phi_1 scaling
        delta_phi_1_small = 0.001
        delta_phi_1_large = 0.002
        
        phase_small = self.generator.phase_ppE_corrections(
            f_test, M_total, delta_phi_1=delta_phi_1_small, delta_phi_1p5=0.0
        )
        phase_large = self.generator.phase_ppE_corrections(
            f_test, M_total, delta_phi_1=delta_phi_1_large, delta_phi_1p5=0.0
        )
        
        # Should scale linearly
        expected_ratio = delta_phi_1_large / delta_phi_1_small
        actual_ratio = phase_large[0] / phase_small[0]
        
        assert abs(actual_ratio - expected_ratio) < 1e-10, "ppE parameters don't scale linearly"
    
    def test_frequency_dependence(self):
        """Test correct frequency dependence of ppE terms."""
        
        # Test frequencies with known ratio
        f1, f2 = 50.0, 100.0  # 2:1 ratio
        M_total = 35.0
        
        # +1PN term should scale as f^(-1)
        phase1_1PN = self.generator.phase_ppE_corrections(
            np.array([f1]), M_total, delta_phi_1=1.0, delta_phi_1p5=0.0
        )[0]
        phase2_1PN = self.generator.phase_ppE_corrections(
            np.array([f2]), M_total, delta_phi_1=1.0, delta_phi_1p5=0.0
        )[0]
        
        expected_ratio_1PN = f2 / f1  # Should be 2.0
        actual_ratio_1PN = phase1_1PN / phase2_1PN
        
        assert abs(actual_ratio_1PN - expected_ratio_1PN) < 1e-6, "1PN frequency dependence incorrect"
        
        # +1.5PN term should scale as f^(-2/3)
        phase1_1p5PN = self.generator.phase_ppE_corrections(
            np.array([f1]), M_total, delta_phi_1=0.0, delta_phi_1p5=1.0
        )[0]
        phase2_1p5PN = self.generator.phase_ppE_corrections(
            np.array([f2]), M_total, delta_phi_1=0.0, delta_phi_1p5=1.0
        )[0]
        
        expected_ratio_1p5PN = (f2 / f1)**(2/3)  # Should be 2^(2/3) ≈ 1.587
        actual_ratio_1p5PN = phase1_1p5PN / phase2_1p5PN
        
        assert abs(actual_ratio_1p5PN - expected_ratio_1p5PN) < 1e-6, "1.5PN frequency dependence incorrect"
    
    def test_mass_dependence(self):
        """Test correct mass dependence of ppE terms."""
        
        f_test = np.array([100.0])
        M1, M2 = 30.0, 60.0  # 2:1 ratio
        
        phase1 = self.generator.phase_ppE_corrections(
            f_test, M1, delta_phi_1=1.0, delta_phi_1p5=0.0
        )[0]
        phase2 = self.generator.phase_ppE_corrections(
            f_test, M2, delta_phi_1=1.0, delta_phi_1p5=0.0
        )[0]
        
        # ppE corrections should scale as M^(-1) for 1PN and M^(-2/3) for 1.5PN
        expected_ratio = M2 / M1  # For 1PN term
        actual_ratio = phase1 / phase2
        
        assert abs(actual_ratio - expected_ratio) < 1e-6, "Mass dependence incorrect"
    
    def test_validation_method(self):
        """Test the built-in validation method."""
        
        validation_passed = self.generator.validate_ppE_implementation()
        
        assert validation_passed, "Built-in validation should pass"


# Integration test
def test_ppE_integration():
    """Integration test for complete ppE pipeline."""
    
    generator = ppEWaveformGenerator()
    
    # Representative parameters
    f_grid = np.logspace(1.3, 2.5, 100)  # 20-300 Hz
    m1, m2 = 30.0, 25.0
    delta_phi_1, delta_phi_1p5 = 0.01, -0.015
    
    # Generate waveform
    h_plus, h_cross = generator.generate_ppE_waveform(
        f_grid, m1, m2, delta_phi_1=delta_phi_1, delta_phi_1p5=delta_phi_1p5
    )
    
    # Basic sanity checks
    assert len(h_plus) == len(f_grid)
    assert len(h_cross) == len(f_grid)
    assert np.all(np.isfinite(h_plus[f_grid > 0]))
    assert np.all(np.isfinite(h_cross[f_grid > 0]))
    
    # Check that signal has reasonable amplitude structure
    # (should decrease with frequency in inspiral)
    mid_idx = len(f_grid) // 2
    assert np.abs(h_plus[mid_idx//2]) > np.abs(h_plus[mid_idx]), "Amplitude should decrease with frequency"


if __name__ == "__main__":
    # Run tests if called directly
    pytest.main([__file__, "-v"])