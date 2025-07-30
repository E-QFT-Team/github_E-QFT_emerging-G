"""
Emergent metric from commutator distance calculation.

This module implements the emergent distance metric d(x_i, x_j)^2 according to Eq. 81:
d(x_i, x_j)^2 = Tr([π_x_i, π_x_j][π_x_i, π_x_j]†)

where π_x are local projectors given by Eq. (69) (Fourier-localization).
"""

import numpy as np
from functools import lru_cache
from typing import Dict, Union
import scipy.sparse as sp

# SIMPLIFIED: Direct calculation for flat space scaling
# No redundant correction factors - _FLAT_PREFAC = 1.0 exactly
_FLAT_PREFAC = 1.0  # Simplified from redundant calculation

# Global dictionary to store projectors
# projectors[index] -> np.ndarray or scipy.sparse.csr_matrix (D x D hermitian and idempotent)
projectors: Dict[int, Union[np.ndarray, sp.csr_matrix]] = {}


@lru_cache(maxsize=None)
def commutator_metric(i: int, j: int) -> float:
    """
    Calculate emergent distance d(x_i, x_j)^2 according to Eq. 81.
    
    The distance is computed as the trace of the commutator squared:
    d(x_i, x_j)^2 = Tr([π_i, π_j][π_i, π_j]†) = ||[π_i, π_j]||²_F
    
    Parameters
    ----------
    i : int
        Index of first site
    j : int
        Index of second site
        
    Returns
    -------
    float
        Squared distance d_ij^2 (non-negative scalar)
        Returns 0.0 if and only if [π_i, π_j] = 0
        
    Raises
    ------
    KeyError
        If projectors for indices i or j are not found
    """
    if i not in projectors:
        raise KeyError(f"Projector for index {i} not found")
    if j not in projectors:
        raise KeyError(f"Projector for index {j} not found")
    
    Pi = projectors[i]
    Pj = projectors[j]
    
    # Compute commutator [Pi, Pj] = Pi @ Pj - Pj @ Pi
    if sp.issparse(Pi) or sp.issparse(Pj):
        # Handle sparse matrices
        C = Pi.dot(Pj) - Pj.dot(Pi)
    else:
        # Handle dense matrices
        C = Pi @ Pj - Pj @ Pi
    
    # Compute ||C||²_F = Tr(C†C) = <C, C>
    if sp.issparse(C):
        # For sparse matrices, convert to dense for vdot calculation
        C_dense = C.toarray()
        dist2 = np.vdot(C_dense, C_dense).real
    else:
        # For dense matrices
        dist2 = np.vdot(C, C).real
    
    # Note: _FLAT_PREFAC = 1.0 for sigma = 1/sqrt(2), so no correction needed
    # dist2 *= _FLAT_PREFAC
    
    return float(dist2)


def clear_cache():
    """Clear the LRU cache for commutator_metric function."""
    commutator_metric.cache_clear()


def set_projectors(new_projectors: Dict[int, Union[np.ndarray, sp.csr_matrix]]):
    """
    Set the global projectors dictionary and clear cache.
    
    Parameters
    ----------
    new_projectors : dict
        Dictionary mapping site indices to projector matrices
        Each projector should be a D x D hermitian and idempotent matrix
    """
    global projectors
    projectors = new_projectors.copy()
    clear_cache()


def validate_projector(P: Union[np.ndarray, sp.csr_matrix], tol: float = 1e-10) -> bool:
    """
    Validate that a matrix is a proper projector (hermitian and idempotent).
    
    Parameters
    ----------
    P : np.ndarray or scipy.sparse.csr_matrix
        Matrix to validate
    tol : float
        Numerical tolerance for validation
        
    Returns
    -------
    bool
        True if P is hermitian and idempotent within tolerance
    """
    if sp.issparse(P):
        P_dense = P.toarray()
    else:
        P_dense = P
    
    # Check hermiticity: P = P†
    is_hermitian = np.allclose(P_dense, P_dense.conj().T, atol=tol)
    
    # Check idempotency: P² = P
    if sp.issparse(P):
        P2 = P.dot(P).toarray()
    else:
        P2 = P_dense @ P_dense
    
    is_idempotent = np.allclose(P2, P_dense, atol=tol)
    
    return is_hermitian and is_idempotent


def get_cache_info():
    """Get cache statistics for the commutator_metric function."""
    return commutator_metric.cache_info()