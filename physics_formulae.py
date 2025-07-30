"""
Symbolic physics formulae for emergent gravity binary corrections.

This module contains the fundamental physics equations in symbolic form,
keeping the physics transparent and separating it from system-specific parameters.
"""

import sympy as sp
import numpy as np

# Define symbolic variables
e, M1, M2, beta = sp.symbols('e M1 M2 beta', positive=True, real=True)
G, c, a, alpha = sp.symbols('G c a alpha', positive=True, real=True)

# Derived quantities
M_asym = abs(M1 - M2) / (M1 + M2)
M_total = M1 + M2

# Binary correction factor (Eq. 1)
C_corr = 1 + beta * M_asym * (2*e - 1)

# Periastron advance (leading GR term)
Delta_phi = 6*sp.pi*G*M_total / (c**2 * a * (1 - e**2))

# Orbital radius as function of true anomaly
nu = sp.symbols('nu', real=True)
r_orbital = a * (1 - e**2) / (1 + e * sp.cos(nu))

# Orbit averaging integral (symbolic)
f_r = sp.Function('f')
orbit_average_symbolic = sp.integrate(f_r(r_orbital), (nu, 0, 2*sp.pi)) / (2*sp.pi)

# Theoretical beta expansion around e = 1/2
e_half = sp.Rational(1, 2)
f_avg = sp.Function('f_avg')
beta_expansion = f_avg(e_half) + sp.diff(f_avg(e), e).subs(e, e_half) * (e - e_half)

# Export lambdified versions for numerical computation
correction_factor_numeric = sp.lambdify((beta, M1, M2, e), C_corr, 'numpy')
mass_asymmetry_numeric = sp.lambdify((M1, M2), M_asym, 'numpy')
periastron_advance_numeric = sp.lambdify((G, c, M1, M2, a, e), Delta_phi, 'numpy')
orbital_radius_numeric = sp.lambdify((a, e, nu), r_orbital, 'numpy')

# Physical constants (SI units)
G_SI = 6.67430e-11  # m³ kg⁻¹ s⁻²
c_SI = 299792458    # m s⁻¹
M_sun = 1.98847e30  # kg

# Single source of truth for beta coefficient
# From theoretical lattice-consistent orbit averaging: β = d⟨f⟩/de|_{e=1/2}
beta_coefficient = -0.033

def get_beta_coefficient():
    """
    Get the current beta coefficient value.
    
    Returns
    -------
    float
        Beta coefficient used in E-QFT corrections
    """
    return beta_coefficient

def set_beta_coefficient(new_beta):
    """
    Set a new beta coefficient value for testing purposes.
    
    Parameters
    ----------
    new_beta : float
        New beta coefficient value
    """
    global beta_coefficient
    beta_coefficient = float(new_beta)
    print(f"🔧 Beta coefficient updated to: {beta_coefficient}")

def correction_factor(M1, M2, e):
    """
    Compute binary correction factor C for gravitational effects.
    
    The correction factor accounts for mass asymmetry and orbital eccentricity:
    C = 1 + β * |M₁-M₂|/(M₁+M₂) * (2e-1)
    
    Parameters
    ----------
    M1 : float
        Mass of primary object (in solar masses or any consistent units)
    M2 : float
        Mass of secondary object (same units as M1)
    e : float
        Orbital eccentricity (0 ≤ e < 1)
        
    Returns
    -------
    float
        Correction factor C
    """
    # Use the numerical version with current beta
    return correction_factor_numeric(beta_coefficient, M1, M2, e)

def display_formulae():
    """
    Display all symbolic formulae for reference.
    """
    print("=== SYMBOLIC PHYSICS FORMULAE ===")
    print(f"Mass asymmetry: M_asym = {M_asym}")
    print(f"Correction factor: C = {C_corr}")
    print(f"Periastron advance: Δφ = {Delta_phi}")
    print(f"Orbital radius: r(ν) = {r_orbital}")
    print(f"Beta expansion: ⟨f⟩(e) ≈ {beta_expansion}")

def validate_formulae():
    """
    Validate formulae with known limits and test cases.
    
    Returns
    -------
    bool
        True if all validation tests pass
    """
    print("=== FORMULA VALIDATION ===")
    
    # Test 1: Equal masses should give M_asym = 0
    test_asym = mass_asymmetry_numeric(1.4, 1.4)
    print(f"Equal masses (1.4, 1.4): M_asym = {test_asym:.10f}")
    assert abs(test_asym) < 1e-10, "Equal masses should give zero asymmetry"
    
    # Test 2: e = 0.5 should make (2e - 1) = 0, so C = 1
    test_corr = correction_factor_numeric(1.0, 1.4, 1.3, 0.5)
    print(f"e = 0.5: C = {test_corr:.10f}")
    assert abs(test_corr - 1.0) < 1e-10, "e = 0.5 should give C = 1"
    
    # Test 3: Circular orbit (e = 0) should give r = a
    test_radius = orbital_radius_numeric(1.0, 0.0, np.pi/4)
    print(f"Circular orbit: r = {test_radius:.10f}")
    assert abs(test_radius - 1.0) < 1e-10, "Circular orbit should give r = a"
    
    # Test 4: Periastron (ν = 0) for eccentric orbit
    test_peri = orbital_radius_numeric(1.0, 0.6, 0.0)
    expected_peri = 1.0 * (1 - 0.6**2) / (1 + 0.6)
    print(f"Periastron (e=0.6): r = {test_peri:.6f}, expected = {expected_peri:.6f}")
    assert abs(test_peri - expected_peri) < 1e-10, "Periastron formula incorrect"
    
    print("✅ All validation tests passed")
    return True

if __name__ == "__main__":
    display_formulae()
    print()
    validate_formulae()