"""
Lattice scale calibration for E-QFT quantum gravity.

Stores the calibration choice for lattice spacing 'a' that converts
from lattice-dimensionless units to SI units.

This is the single source of truth for the lattice scale calibration.
"""

# Lattice spacing calibration (metres)
# This value converts from lattice units to SI units
# See manuscript §3 for derivation and justification
a_SI = 9.23e-35  # metres

# Lattice calibration metadata
calibration_info = {
    "value": a_SI,
    "units": "metres",
    "source": "manuscript §3",
    "description": "Lattice spacing that converts lattice-dimensionless G_eff to SI units",
    "uncertainty": "TBD - requires lattice systematics analysis"
}

def get_lattice_spacing():
    """
    Get the calibrated lattice spacing.
    
    Returns
    -------
    float
        Lattice spacing in metres
    """
    return a_SI

def display_calibration():
    """Display lattice calibration information."""
    print("📏 LATTICE SCALE CALIBRATION")
    print("=" * 30)
    print(f"Lattice spacing: a = {a_SI:.2e} m")
    print(f"Source: {calibration_info['source']}")
    print(f"Description: {calibration_info['description']}")

if __name__ == "__main__":
    display_calibration()