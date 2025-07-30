"""
Unit tests for E-QFT integration into ppE framework.

Tests the lattice → SI unit bridge implementation as specified in the task.
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from constants import G_N_SI
from geff_conversion import Geff_SI, validate_geff_units
from physics_formulae import correction_factor_numeric
from physics_formulae import get_beta_coefficient
from ppE_waveform_module import ppEWaveformGenerator


class TestGeffUnits:
    """Test effective gravitational constant units and conversion."""
    
    def test_geff_units(self):
        """Test that G_eff is in the correct physical range."""
        G = Geff_SI()
        # Must straddle Newton's G as specified in task
        assert 6e-11 < G < 7e-11, f"G_eff = {G:.2e} outside reasonable range"
    
    def test_geff_validation(self):
        """Test G_eff validation function."""
        assert validate_geff_units() == True
    
    def test_geff_relative_to_newton(self):
        """Test G_eff relative to Newton's constant."""
        G_eff = Geff_SI()
        ratio = G_eff / G_N_SI
        # Should be within factor of 2 of Newton's G
        assert 0.5 < ratio < 2.0, f"G_eff/G_N = {ratio:.3f} too far from unity"


class TestCorrectionFactor:
    """Test binary correction factor implementation."""
    
    def test_correction_factor_b1913(self):
        """Test correction factor for PSR B1913+16."""
        beta = get_beta_coefficient()  # Should be -0.033
        C = correction_factor_numeric(beta, 1.441, 1.387, 0.617155)
        
        # With beta = -0.033, expect correction near 1.0
        assert 0.999 < C < 1.001, f"Correction factor C = {C:.6f} outside expected range"
    
    def test_correction_factor_circular_equal(self):
        """Test correction factor for circular equal-mass system."""
        beta = get_beta_coefficient()
        C = correction_factor_numeric(beta, 1.4, 1.4, 0.0)
        
        # Equal masses, circular orbit should give C = 1
        assert abs(C - 1.0) < 1e-10, f"Equal masses should give C = 1, got {C:.10f}"
    
    def test_correction_factor_e_half(self):
        """Test correction factor at e = 0.5."""
        beta = get_beta_coefficient()
        C = correction_factor_numeric(beta, 1.4, 1.3, 0.5)
        
        # e = 0.5 should give (2e-1) = 0, so C = 1
        assert abs(C - 1.0) < 1e-10, f"e = 0.5 should give C = 1, got {C:.10f}"
    
    def test_beta_coefficient_range(self):
        """Test that beta coefficient is in expected range."""
        beta = get_beta_coefficient()
        
        # Our theoretical derivation gives beta = -0.033
        assert -0.1 < beta < 0.1, f"Beta = {beta} outside reasonable range"


class TestppEIntegration:
    """Test E-QFT integration into ppE waveform framework."""
    
    def test_generator_initialization(self):
        """Test that ppE generator initializes with E-QFT corrections."""
        gen = ppEWaveformGenerator()
        
        # Should have both Newton's G and E-QFT G_eff
        assert hasattr(gen, 'G_newton')
        assert hasattr(gen, 'G_eff')
        assert hasattr(gen, 'beta')
        
        # G_eff should be different from Newton's G
        assert gen.G_eff != gen.G_newton
    
    def test_effective_gravity_computation(self):
        """Test effective gravity computation."""
        gen = ppEWaveformGenerator()
        
        # Test with PSR B1913+16 parameters
        G_corrected = gen.compute_effective_gravity(1.441, 1.387, 0.617155)
        
        # Should be close to G_eff since correction is small
        assert abs(G_corrected - gen.G_eff) / gen.G_eff < 0.01
    
    def test_system_card_loading(self):
        """Test loading system parameters from card."""
        gen = ppEWaveformGenerator()
        
        # Test with PSR B1913+16 card
        params = gen.load_system_card('cards/PSR_B1913+16.yaml')
        
        # Check required parameters are present and numeric
        assert 'M1' in params and isinstance(params['M1'], (int, float))
        assert 'M2' in params and isinstance(params['M2'], (int, float))
        assert 'e' in params and isinstance(params['e'], (int, float))
    
    @pytest.mark.slow
    def test_eqft_waveform_generation(self):
        """Test E-QFT waveform generation."""
        gen = ppEWaveformGenerator()
        test_frequencies = np.logspace(1, 3, 100)  # 10 Hz to 1 kHz
        
        try:
            # Generate waveform with E-QFT corrections
            h_plus, h_cross = gen.frequency_domain_waveform_eqft(
                'cards/PSR_B1913+16.yaml',
                test_frequencies,
                delta_phi_1=0.01,
                delta_phi_1p5=0.02
            )
            
            # Check output format
            assert len(h_plus) == len(test_frequencies)
            assert len(h_cross) == len(test_frequencies)
            assert np.all(np.isfinite(h_plus))
            assert np.all(np.isfinite(h_cross))
            
        except FileNotFoundError:
            pytest.skip("Parameter card not found - test requires cards/PSR_B1913+16.yaml")


class TestTaskSpecification:
    """Test compliance with task specification."""
    
    def test_constants_centralization(self):
        """Test that physical constants are centralized."""
        from constants import c_SI, G_N_SI, hbar_SI, m_P_SI
        
        # Constants should be properly defined
        assert c_SI > 0
        assert G_N_SI > 0
        assert hbar_SI > 0
        assert m_P_SI > 0
    
    def test_lattice_scale_single_source(self):
        """Test that lattice scale is in single location."""
        from lattice_scale import a_SI, get_lattice_spacing
        
        # Should have single source of truth for lattice spacing
        assert a_SI > 0
        assert get_lattice_spacing() == a_SI
    
    def test_geff_conversion_formula(self):
        """Test that G_eff conversion follows Eq.(1) from task."""
        from constants import c_SI, m_P_SI
        from lattice_scale import a_SI
        from geff_conversion import G_eff_lat, Geff_SI
        
        # Manually compute Eq.(1): G_eff,SI = G_eff,lat × a × c² / m_P
        expected = G_eff_lat * a_SI * c_SI**2 / m_P_SI
        actual = Geff_SI()
        
        # Should match exactly
        assert abs(actual - expected) / expected < 1e-10
    
    def test_binary_correction_integration(self):
        """Test that binary correction is properly integrated."""
        gen = ppEWaveformGenerator()
        
        # Should use theoretical beta from gravity_corrections
        from physics_formulae import get_beta_coefficient
        assert gen.beta == get_beta_coefficient()
        
        # Effective gravity should include binary correction
        G1 = gen.compute_effective_gravity(1.4, 1.4, 0.0)  # No correction
        G2 = gen.compute_effective_gravity(2.0, 1.0, 0.8)  # Large correction
        
        # Should be different due to binary correction
        assert G1 != G2


def run_integration_tests():
    """Run all E-QFT integration tests."""
    print("🧪 E-QFT INTEGRATION TESTS")
    print("=" * 40)
    
    test_classes = [
        TestGeffUnits,
        TestCorrectionFactor,
        TestppEIntegration,
        TestTaskSpecification
    ]
    
    total_tests = 0
    passed_tests = 0
    
    for test_class in test_classes:
        print(f"\n📋 {test_class.__name__}")
        instance = test_class()
        
        # Special setup for ppE tests
        if hasattr(instance, 'setUp'):
            instance.setUp()
        
        for method_name in dir(instance):
            if method_name.startswith('test_'):
                total_tests += 1
                try:
                    method = getattr(instance, method_name)
                    method()
                    print(f"  ✅ {method_name}")
                    passed_tests += 1
                except Exception as e:
                    print(f"  ❌ {method_name}: {e}")
    
    print(f"\n🎯 RESULTS: {passed_tests}/{total_tests} tests passed")
    
    if passed_tests == total_tests:
        print("🎉 ALL E-QFT INTEGRATION TESTS PASSED!")
        return True
    else:
        print("⚠️  Some integration tests failed")
        return False


if __name__ == "__main__":
    success = run_integration_tests()
    exit(0 if success else 1)