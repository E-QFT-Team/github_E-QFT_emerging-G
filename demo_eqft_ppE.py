#!/usr/bin/env python3
"""
Demonstration of E-QFT integration into ppE framework.

Shows how the lattice → SI unit bridge brings theoretical β coefficient
into gravitational waveform generation with proper physics separation.
"""

import numpy as np
import matplotlib.pyplot as plt
from ppE_waveform_module import ppEWaveformGenerator

def demo_eqft_integration():
    """Demonstrate complete E-QFT integration."""
    print("🚀 E-QFT ppE INTEGRATION DEMONSTRATION")
    print("=" * 50)
    
    # Initialize generator with E-QFT corrections
    gen = ppEWaveformGenerator()
    
    print(f"\n📊 THEORETICAL FOUNDATION")
    print(f"Lattice β coefficient: {gen.beta:.6f}")
    print(f"Newton's G: {gen.G_newton:.2e} m³/(kg·s²)")
    print(f"E-QFT G_eff: {gen.G_eff:.2e} m³/(kg·s²)")
    print(f"Ratio G_eff/G_N: {gen.G_eff/gen.G_newton:.6f}")
    
    # Test with different binary systems
    systems = ["PSR_B1913+16", "Mercury"]
    
    for system in systems:
        print(f"\n🔬 SYSTEM: {system}")
        print("-" * 30)
        
        try:
            # Load system parameters
            params = gen.load_system_card(f"cards/{system}.yaml")
            M1, M2, e = params['M1'], params['M2'], params['e']
            
            # Compute effective gravity
            G_corrected = gen.compute_effective_gravity(M1, M2, e)
            correction_factor = G_corrected / gen.G_eff
            
            print(f"Masses: M₁={M1:.3f} M☉, M₂={M2:.3f} M☉")
            print(f"Eccentricity: e={e:.6f}")
            print(f"Binary correction: C={correction_factor:.6f}")
            print(f"Corrected G: {G_corrected:.6e} m³/(kg·s²)")
            
            # Show ppE parameter effects
            print(f"\n📈 ppE PARAMETER EFFECTS")
            
            # Test frequency range (LIGO band)
            f_array = np.logspace(1, 3, 1000)  # 10 Hz to 1 kHz
            
            # Generate waveforms with different ppE parameters
            ppE_cases = [
                {"delta_phi_1": 0.0, "delta_phi_1p5": 0.0, "label": "GR + E-QFT"},
                {"delta_phi_1": 0.01, "delta_phi_1p5": 0.0, "label": "GR + E-QFT + δφ₁"},
                {"delta_phi_1": 0.0, "delta_phi_1p5": 0.02, "label": "GR + E-QFT + δφ₁.₅"},
                {"delta_phi_1": 0.01, "delta_phi_1p5": 0.02, "label": "GR + E-QFT + δφ₁ + δφ₁.₅"}
            ]
            
            waveforms = []
            for case in ppE_cases:
                try:
                    h_plus, h_cross = gen.frequency_domain_waveform_eqft(
                        f"cards/{system}.yaml",
                        f_array,
                        delta_phi_1=case["delta_phi_1"],
                        delta_phi_1p5=case["delta_phi_1p5"]
                    )
                    waveforms.append({
                        "f": f_array,
                        "h_plus": h_plus,
                        "h_cross": h_cross,
                        "label": case["label"]
                    })
                    print(f"  ✅ Generated: {case['label']}")
                except Exception as e:
                    print(f"  ❌ Failed: {case['label']} - {e}")
            
            # Create comparison plot
            if waveforms:
                create_waveform_comparison_plot(waveforms, system)
                
        except FileNotFoundError:
            print(f"❌ Card not found: cards/{system}.yaml")
        except Exception as e:
            print(f"❌ Error processing {system}: {e}")
    
    print(f"\n✅ DEMONSTRATION COMPLETE!")
    print(f"Key achievements:")
    print(f"• Lattice β = {gen.beta:.6f} integrated into ppE framework")
    print(f"• G_eff from theoretical lattice → SI conversion")
    print(f"• Binary corrections from parameter cards")
    print(f"• Full waveform generation with E-QFT + ppE")

def create_waveform_comparison_plot(waveforms, system_name):
    """Create comparison plot of waveforms with different corrections."""
    plt.figure(figsize=(12, 8))
    
    # Plot amplitude comparison
    plt.subplot(2, 1, 1)
    for wf in waveforms:
        amplitude = np.abs(wf["h_plus"] + 1j * wf["h_cross"])
        plt.loglog(wf["f"], amplitude, label=wf["label"], linewidth=2)
    
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("|h(f)| [strain]")
    plt.title(f"{system_name}: Amplitude Comparison (E-QFT + ppE)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot phase difference relative to GR+E-QFT
    plt.subplot(2, 1, 2)
    if len(waveforms) > 1:
        reference_phase = np.angle(waveforms[0]["h_plus"] + 1j * waveforms[0]["h_cross"])
        
        for i, wf in enumerate(waveforms[1:], 1):
            phase = np.angle(wf["h_plus"] + 1j * wf["h_cross"])
            phase_diff = phase - reference_phase
            
            # Unwrap phase differences
            phase_diff = np.unwrap(phase_diff)
            
            plt.semilogx(wf["f"], phase_diff, label=wf["label"], linewidth=2)
    
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Δφ [rad] (relative to GR+E-QFT)")
    plt.title(f"{system_name}: Phase Deviations from ppE Parameters")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_filename = f"eqft_ppE_comparison_{system_name.lower()}.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"  📊 Plot saved: {plot_filename}")
    plt.close()

def validate_physics_consistency():
    """Validate that physics is consistent across the integration."""
    print(f"\n🔬 PHYSICS VALIDATION")
    print("-" * 25)
    
    gen = ppEWaveformGenerator()
    
    # Test 1: e=0.5 should give C=1
    G_half = gen.compute_effective_gravity(1.4, 1.3, 0.5)
    ratio_half = G_half / gen.G_eff
    print(f"e=0.5 test: C = {ratio_half:.10f} (should be 1.0)")
    assert abs(ratio_half - 1.0) < 1e-10, "e=0.5 physics violation!"
    
    # Test 2: Equal masses should give C=1
    G_equal = gen.compute_effective_gravity(1.4, 1.4, 0.617)
    ratio_equal = G_equal / gen.G_eff
    print(f"Equal masses: C = {ratio_equal:.10f} (should be 1.0)")
    assert abs(ratio_equal - 1.0) < 1e-10, "Equal mass physics violation!"
    
    # Test 3: G_eff should be close to G_Newton
    ratio_geff = gen.G_eff / gen.G_newton
    print(f"G_eff/G_N = {ratio_geff:.6f} (should be ≈ 1)")
    assert 0.5 < ratio_geff < 2.0, "G_eff unreasonable!"
    
    print("✅ All physics consistency tests passed!")

if __name__ == "__main__":
    # Run demonstration
    demo_eqft_integration()
    
    # Validate physics
    validate_physics_consistency()
    
    print(f"\n🎉 E-QFT ppE INTEGRATION SUCCESSFUL!")
    print(f"The theoretical β = -0.033 is now properly integrated")
    print(f"into the production ppE waveform framework!")