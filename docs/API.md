# API Reference

## Core Modules

### projector_factory

#### `generate_fourier_projectors(lattice_shape, localization_width=1.0, lambda_param=1.0, use_sparse=False, dimension=64)`

Generate local projectors using Fourier-localization.

**Parameters:**
- `lattice_shape` (tuple): Shape of the lattice, e.g., (8, 8, 8) for 8³
- `localization_width` (float): Gaussian localization parameter σ
- `lambda_param` (float): Coupling parameter λ
- `use_sparse` (bool): Return sparse matrices if True
- `dimension` (int): Matrix dimension D (D×D projectors)

**Returns:**
- `dict`: Mapping from site indices to projector matrices

**Example:**
```python
from src.projector_factory import generate_fourier_projectors

projectors = generate_fourier_projectors(
    lattice_shape=(8, 8, 8),
    localization_width=1.0,
    dimension=64
)
```

#### `generate_1d_test_projectors(num_sites, lambda_param=0.01, dimension=32, use_sparse=False)`

Generate projectors for 1D flat space test.

**Parameters:**
- `num_sites` (int): Number of 1D lattice sites
- `lambda_param` (float): Coupling parameter (kept for compatibility)
- `dimension` (int): Matrix dimension
- `use_sparse` (bool): Return sparse matrices if True

**Returns:**
- `dict`: Projectors that yield exact flat space scaling d²(i,j) = (i-j)²

### metric

#### `set_projectors(projector_dict)`

Set the global projector dictionary for metric calculations.

**Parameters:**
- `projector_dict` (dict): Mapping from site indices to projector matrices

#### `commutator_metric(i, j)`

Compute commutator metric d²(i,j) = Tr([Πᵢ, Πⱼ]†[Πᵢ, Πⱼ]).

**Parameters:**
- `i` (int): First site index
- `j` (int): Second site index

**Returns:**
- `float`: Metric distance squared

**Example:**
```python
from src.metric import set_projectors, commutator_metric

set_projectors(projectors)
distance = commutator_metric(0, 1)
```

#### `clear_cache()`

Clear the internal cache for metric calculations.

#### `validate_projector(P, hermitian_tol=1e-12)`

Validate that a matrix is a valid projector (Hermitian).

**Parameters:**
- `P` (np.ndarray): Matrix to validate
- `hermitian_tol` (float): Tolerance for Hermiticity check

**Returns:**
- `bool`: True if valid projector

### entropy_tools

#### `von_neumann_entropy(density_matrix, eigenvalue_threshold=1e-12)`

Compute von Neumann entropy S = -Tr(ρ log ρ).

**Parameters:**
- `density_matrix` (np.ndarray): Density matrix ρ
- `eigenvalue_threshold` (float): Cutoff for zero eigenvalues

**Returns:**
- `float`: Von Neumann entropy

**Example:**
```python
from src.entropy_tools import von_neumann_entropy

# Create density matrix
rho = np.array([[0.7, 0], [0, 0.3]])
entropy = von_neumann_entropy(rho)
```

## Analysis Modules

### final_comprehensive_validation

#### `create_final_validation_summary()`

Generate comprehensive validation summary with all results.

**Returns:**
- `dict`: Complete validation results including edge profile, radius convergence, Protocol B, and normalization

#### `create_final_validation_plots()`

Create final comprehensive validation plots.

### protocol_b_stabilization

#### `test_mu_linearity()`

Test κ ∝ μ linearity as specified in Protocol B.

**Returns:**
- `dict`: Linearity test results including slope, intercept, and pass/fail status

#### `fit_kappa(mu, lattice_size=10)`

Fit κ parameter for given mass-defect μ.

**Parameters:**
- `mu` (float): Mass-defect parameter
- `lattice_size` (int): Lattice size

**Returns:**
- `dict`: Fit results including κ value, error, and R²

### metric_normalization

#### `NormalizedMetric.normalized_commutator_metric(i, j)`

Compute normalized metric with d̄²_nn = 1.

**Parameters:**
- `i` (int): First site index  
- `j` (int): Second site index

**Returns:**
- `float`: Normalized metric distance

### geff_rc_scan

#### `geff_for_radius(lattice_n=8, max_rc=4, block_ratio=0.5)`

Compute G_eff for extended connectivity up to max_rc.

**Parameters:**
- `lattice_n` (int): Lattice size (n³)
- `max_rc` (int): Maximum connectivity radius
- `block_ratio` (float): Block size ratio for region definition

**Returns:**
- `dict`: G_eff results with connectivity analysis

## Utility Functions

### Graph Construction

Most analysis modules use NetworkX for graph construction:

```python
import networkx as nx

G = nx.Graph()
for i in range(total_sites):
    for j in range(i+1, total_sites):
        if chebyshev_distance(i, j) <= max_radius:
            weight = commutator_metric(i, j)
            G.add_edge(i, j, weight=weight)
```

### Distance Calculations

Chebyshev distance for cubic lattices:

```python
def chebyshev_distance(i, j, lattice_shape):
    coords_i = np.unravel_index(i, lattice_shape)
    coords_j = np.unravel_index(j, lattice_shape)
    return max(abs(coords_i[k] - coords_j[k]) for k in range(len(lattice_shape)))
```

### Common Patterns

#### Basic Workflow

```python
# 1. Generate projectors
projectors = generate_fourier_projectors((8, 8, 8), dimension=64)

# 2. Set for metric calculations  
set_projectors(projectors)

# 3. Compute distances
d_squared = commutator_metric(i, j)

# 4. Build graph
G = nx.Graph()
# ... add edges with weights

# 5. Apply min-cut algorithm
cut_value, partition = nx.minimum_cut(G, source_nodes, sink_nodes)
```

#### Error Handling

All functions include appropriate error handling:

```python
try:
    result = commutator_metric(i, j)
except ValueError as e:
    print(f"Invalid site indices: {e}")
```

#### Memory Management

For large lattices, use sparse matrices:

```python
projectors = generate_fourier_projectors(
    (16, 16, 16), 
    use_sparse=True,  # Reduce memory usage
    dimension=64
)
```