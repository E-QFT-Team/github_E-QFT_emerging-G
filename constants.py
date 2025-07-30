"""
Physical constants for E-QFT quantum gravity calculations.

Uses scipy.constants for well-known physical constants to avoid hard-coding values.
"""

import numpy as np
import scipy.constants as const

# Fundamental SI constants from scipy
c_SI = const.c               # m s⁻¹ (speed of light)
G_N_SI = const.G             # m³ kg⁻¹ s⁻² (Newton's gravitational constant)
hbar_SI = const.hbar         # J s (reduced Planck constant)

# Derived constants
m_P_SI = np.sqrt(hbar_SI * c_SI / G_N_SI)  # kg (Planck mass)

# Astronomical constants
M_sun_SI = 1.98847e30        # kg (solar mass - not in scipy.constants)

# Display constants for verification
def display_constants():
    """Display all physical constants for verification."""
    print("🔬 FUNDAMENTAL PHYSICAL CONSTANTS")
    print("=" * 40)
    print(f"Speed of light:    c = {c_SI:.8e} m/s")
    print(f"Newton's G:        G_N = {G_N_SI:.8e} m³/(kg·s²)")
    print(f"Reduced Planck ℏ:  ℏ = {hbar_SI:.9e} J·s")
    print(f"Planck mass:       m_P = {m_P_SI:.6e} kg")
    print(f"Solar mass:        M_☉ = {M_sun_SI:.5e} kg")

if __name__ == "__main__":
    display_constants()