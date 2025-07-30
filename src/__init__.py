"""
Emergent Gravity: Core source code package.

This package contains the fundamental components for emergent gravity calculations:
- projector_factory: Deterministic rank-1 projector generation
- metric: Quantum commutator metric computation  
- entropy_tools: Von Neumann entropy calculation

Key classes and functions are imported at package level for convenience.
"""

from .projector_factory import (
    generate_fourier_projectors,
    generate_1d_test_projectors,
    create_test_lattice_8cubed,
    create_test_lattice_16cubed,
)

from .metric import (
    commutator_metric,
    set_projectors,
    clear_cache,
    validate_projector,
)

from .entropy_tools import (
    von_neumann_entropy,
)

__version__ = "1.0.0"
__author__ = "Emergent Gravity Research Team"

__all__ = [
    # Projector factory
    "generate_fourier_projectors",
    "generate_1d_test_projectors", 
    "create_test_lattice_8cubed",
    "create_test_lattice_16cubed",
    
    # Metric computation
    "commutator_metric",
    "set_projectors", 
    "clear_cache",
    "validate_projector",
    
    # Entropy calculation
    "von_neumann_entropy",
]