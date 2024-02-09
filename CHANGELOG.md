# Changelog

This file highlights the major changes made to `poissbox`, it is based on the recommendations made
by [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the versioning system follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Implement serial Poisson solver based on compact schemes
- Implement compact Laplacian operator based on div(grad(f))
- Implement compact schemes for gradient and divergence operations
- Implement tridiagonal solver
- Implement matrix-free version of linear system
- Solve linear system and check error
- Assembly of 2nd-order Laplacian operator and evaluation via matrix-vector product
- Solution vector value initialisation (for computing reference RHS, etc.)
- Build a linear system (matrix, solution + RHS vectors) from the mesh object
- Build a 3-D distributed structured mesh using PETSc DMDA
- Add MPI + PETSc integration

### Changed
### Deprecated
### Removed
### Fixed
