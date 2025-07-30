"""
Effective gravitational constant conversion from lattice to SI units.

Implements the lattice → SI unit bridge:
G_eff,SI = G_eff,lat × a × c² / m_P

This provides the absolute scale of gravity for all waveform calculations.
"""

from constants import c_SI, m_P_SI
from lattice_scale import a_SI

# Lattice-dimensionless effective gravitational constant
# From lattice QFT calculations
G_eff_lat = 0.171772

def Geff_SI():
    """
    Convert lattice effective G to SI units.
    
    Implements Eq.(1): G_eff,SI = G_eff,lat × a × c² / m_P
    
    Returns
    -------
    float
        Effective gravitational constant in SI units (m³ kg⁻¹ s⁻²)
    """
    return G_eff_lat * a_SI * c_SI**2 / m_P_SI

def get_lattice_geff():
    """
    Get the lattice-dimensionless G_eff value.
    
    Returns
    -------
    float
        Lattice G_eff (dimensionless)
    """
    return G_eff_lat

def display_conversion():
    """Display G_eff conversion details."""
    G_SI = Geff_SI()
    
    print("⚖️  EFFECTIVE GRAVITATIONAL CONSTANT CONVERSION")
    print("=" * 50)
    print(f"Lattice G_eff:     G_eff,lat = {G_eff_lat}")
    print(f"Lattice spacing:   a = {a_SI:.2e} m")
    print(f"Speed of light:    c = {c_SI:.2e} m/s")
    print(f"Planck mass:       m_P = {m_P_SI:.2e} kg")
    print()
    print(f"SI G_eff:          G_eff,SI = {G_SI:.6e} m³/(kg·s²)")
    print(f"Newton's G:        G_N = 6.67430e-11 m³/(kg·s²)")
    print(f"Ratio G_eff/G_N:   {G_SI/6.67430e-11:.6f}")

def validate_geff_units():
    """
    Validate that G_eff is in the correct range.
    
    Returns
    -------
    bool
        True if G_eff is physically reasonable
    """
    G = Geff_SI()
    G_N = 6.67430e-11
    
    # Should be close to Newton's G but allow for quantum corrections
    valid = 0.5 * G_N < G < 2.0 * G_N
    
    if not valid:
        print(f"⚠️  WARNING: G_eff = {G:.2e} is outside reasonable range")
        print(f"   Expected: {0.5*G_N:.2e} < G_eff < {2.0*G_N:.2e}")
    
    return valid

if __name__ == "__main__":
    display_conversion()
    print()
    print("🧪 VALIDATION")
    print("-" * 20)
    valid = validate_geff_units()
    print(f"G_eff units valid: {'✅' if valid else '❌'}")