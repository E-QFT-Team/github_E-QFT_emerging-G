# Contributing to Emergent Gravity

Thank you for your interest in contributing to this emergent gravity research project! This document provides guidelines for contributing to the codebase.

## Development Philosophy

This project implements theoretical physics calculations that must be:
- **Mathematically rigorous**: All algorithms implement published theoretical results
- **Numerically stable**: Calculations must be robust against numerical precision issues  
- **Scientifically reproducible**: Results must be verifiable and well-documented
- **Computationally efficient**: Code must scale to larger lattice sizes

## How to Contribute

### 1. Types of Contributions

We welcome contributions in several areas:

**Scientific Improvements:**
- Enhanced projector constructions preserving locality
- More efficient metric computation algorithms
- Better finite-size extrapolation methods
- Advanced entropy calculation techniques

**Code Quality:**
- Performance optimizations
- Memory usage improvements
- Better error handling
- Enhanced documentation

**Validation and Testing:**
- Additional unit tests
- Benchmarking against known results
- Cross-validation with other methods
- Edge case testing

**Analysis Tools:**
- Visualization improvements
- Statistical analysis enhancements
- New diagnostic tools
- Better data export formats

### 2. Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/your-username/emergent-gravity.git
   cd emergent-gravity
   ```
3. **Install in development mode**:
   ```bash
   pip install -e .
   pip install -r requirements.txt
   ```
4. **Run tests** to ensure everything works:
   ```bash
   python -m pytest tests/
   ```

### 3. Development Workflow

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the guidelines below

3. **Add tests** for new functionality:
   ```bash
   # Add tests to tests/ directory
   python -m pytest tests/test_your_feature.py
   ```

4. **Run the full test suite**:
   ```bash
   python -m pytest tests/
   ```

5. **Update documentation** if needed

6. **Commit with clear messages**:
   ```bash
   git commit -m "feat: add improved projector locality preservation"
   ```

7. **Push and create a pull request**

## Coding Standards

### Python Style

- Follow **PEP 8** style guidelines
- Use **type hints** for function parameters and returns
- Write **docstrings** for all public functions using NumPy format
- Keep **line length ≤ 88 characters** (Black formatter standard)

Example:
```python
def commutator_metric(i: int, j: int) -> float:
    """
    Compute commutator metric d²(i,j) = Tr([Πᵢ, Πⱼ]†[Πᵢ, Πⱼ]).
    
    Parameters
    ----------
    i : int
        First site index
    j : int
        Second site index
        
    Returns
    -------
    float
        Metric distance squared
        
    Raises
    ------
    ValueError
        If site indices are invalid
    """
```

### Scientific Code Guidelines

1. **No Magic Numbers**: Use named constants for all physical parameters
   ```python
   # Good
   LATTICE_SPACING = 1e-15  # meters
   SPEED_OF_LIGHT = 299792458  # m/s
   
   # Bad
   result = value * 1e-15 * 299792458**2
   ```

2. **Numerical Stability**: Add appropriate tolerances and checks
   ```python
   # Check for numerical issues
   if np.linalg.norm(matrix - matrix.conj().T) > 1e-12:
       raise ValueError("Matrix not Hermitian within tolerance")
   ```

3. **Physical Units**: Always document units and provide conversion utilities
   ```python
   def convert_to_si(G_eff: float, lattice_spacing: float) -> float:
       """Convert G_eff from lattice units to SI units (m³/kg/s²)."""
   ```

4. **Validation**: Include physics-based validation tests
   ```python
   # Test that projectors preserve locality
   def test_locality_preservation():
       # Verify nearby sites have stronger coupling
       assert commutator_metric(i, i+1) < commutator_metric(i, i+10)
   ```

### Performance Guidelines

1. **Use NumPy efficiently**: Prefer vectorized operations over loops
2. **Memory management**: Use sparse matrices for large systems
3. **Caching**: Cache expensive calculations when appropriate
4. **Profiling**: Profile performance-critical sections

## Testing Requirements

### Unit Tests

All new functions must include unit tests:

```python
def test_metric_symmetry():
    """Test that d²(i,j) = d²(j,i)."""
    for i in range(5):
        for j in range(5):
            assert abs(commutator_metric(i, j) - commutator_metric(j, i)) < 1e-12
```

### Physics Validation Tests

Include tests that verify physical correctness:

```python
def test_flat_space_limit():
    """Test exact flat space scaling d²(i,j) = (i-j)²."""
    projectors = generate_1d_test_projectors(10)
    set_projectors(projectors)
    
    for i in range(10):
        for j in range(10):
            expected = (i - j) ** 2
            actual = commutator_metric(i, j)
            assert abs(actual - expected) < 1e-10
```

### Integration Tests

Test complete workflows:

```python
def test_full_protocol_a_workflow():
    """Test complete Protocol A from projectors to G_eff."""
    # Test end-to-end calculation
    # Verify results are physically reasonable
```

## Documentation

### Code Documentation

- **Docstrings**: All public functions need comprehensive docstrings
- **Inline comments**: Explain complex physics or numerical algorithms  
- **Type hints**: Use throughout for better IDE support

### Theoretical Documentation

When adding new theoretical elements:

1. **Update THEORY.md** with mathematical background
2. **Add references** to relevant papers
3. **Include derivations** for non-trivial formulas
4. **Explain physical interpretation**

### Examples

Add working examples for new features:

```python
# examples/new_feature_demo.py
"""
Demonstration of new feature X for emergent gravity calculations.
"""
```

## Scientific Validation

### Required Checks

Before submitting, ensure your changes pass:

1. **Flat space test**: d²(i,j) = (i-j)² exactly
2. **Locality preservation**: Link weights decay with distance
3. **Numerical stability**: Results stable under parameter changes
4. **Physical units**: Dimensional analysis checks out

### Cross-Validation

When possible, validate against:
- Known analytical results
- Previous implementations  
- Independent numerical methods
- Published benchmarks

## Submission Guidelines

### Pull Request Process

1. **Clear title**: Describe the change succinctly
2. **Detailed description**: Explain motivation and approach
3. **Test results**: Include test output and validation results
4. **Performance impact**: Note any performance changes
5. **Breaking changes**: Highlight any API changes

### Review Criteria

Pull requests are evaluated on:

- **Scientific correctness**: Does the physics make sense?
- **Code quality**: Is it well-written and maintainable?
- **Test coverage**: Are all code paths tested?
- **Documentation**: Is it properly documented?
- **Performance**: Does it maintain or improve efficiency?

### Example PR Template

```markdown
## Summary
Brief description of the changes

## Motivation
Why is this change needed?

## Changes Made
- List of specific changes
- New features added
- Bug fixes

## Testing
- Unit tests added/updated
- Physics validation results
- Performance benchmarks

## Documentation
- Updated docstrings
- Theory documentation changes
- New examples added
```

## Questions and Support

- **Scientific questions**: Open an issue with the "question" label
- **Bug reports**: Use the bug report template
- **Feature requests**: Use the feature request template
- **General discussion**: Use GitHub Discussions

## Code of Conduct

This project follows a standard code of conduct:

- Be respectful and inclusive
- Focus on scientific merit
- Provide constructive feedback
- Help maintain a collaborative environment

Thank you for contributing to advancing our understanding of emergent gravity!