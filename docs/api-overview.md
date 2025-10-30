# AlmondFEM API overview

This document captures the solver components that ship with AlmondFEM and the supporting
utilities bundled in the repository. Use it to locate headers, understand each abstraction,
and discover extension points for meshing, boundary conditions, and linear algebra.

## Core headers
The AlmondFEM headers live in `include/almond_fem/` and are delivered through the
`AlmondFEM` CMake interface target. Include only the pieces you need—everything is header
only.

### `mesh.hpp`
- Defines the geometric primitives: `Node`, `Element`, and the owning `Mesh` container.
- `Mesh` validates connectivity, exposes spans for iteration, and offers helper functions
  such as `add_node`, `add_element`, and `node_coordinates`.
- Mesh data can be constructed programmatically or imported from external formats before
  being passed to the solver.

### `problem.hpp`
- Describes the PDE being solved via `ProblemDefinition`.
- Supports concentrated loads (`PointLoad`), uniform source terms, and
  `DirichletBoundary` conditions.
- Returns a `SolveResult` containing the solution vector, residual norm, and iteration count
  when paired with `solve` from `solver.hpp`.

### `solver.hpp`
- Exposes `SolverType` (direct vs. conjugate gradient) and `PreconditionerType` (None,
  Jacobi, IC(0)).
- `SolverOptions` configures tolerances, iteration caps, pivot thresholds, SELL-C-sigma
  slicing, and diagnostic verbosity.
- Implements the element assembly, sparse promotion (COO → CSR), optional SELL/ELL views,
  and solver orchestration.

```cpp
almond::fem::SolverOptions options{};
options.solver = almond::fem::SolverType::ConjugateGradient;
options.preconditioner = almond::fem::PreconditionerType::Jacobi;
options.tolerance = 1e-9;
options.max_iterations = 256;
options.verbose = true;

const auto result = almond::fem::solve(mesh, problem, options);
safe_io::print("Iterations: {} residual: {:.3e}", result.iterations, result.residual_norm);
```

## Assembly workflow
1. Construct a `Mesh` by adding nodes (coordinates) followed by triangular elements referencing
   those nodes. Optionally set per-element conductivity.
2. Populate `ProblemDefinition` with source terms and boundary conditions.
3. Choose solver and preconditioner modes in `SolverOptions`.
4. Call `almond::fem::solve` to receive the nodal solution vector and diagnostics.

## Extension points
- **Custom importers** – wrap mesh creation in helper functions that read your preferred
  file formats; `Mesh` spans simplify zero-copy adapters.
- **New preconditioners** – extend the internal preconditioner helpers inside `solver.hpp`
  and thread the option through `SolverOptions`.
- **Post-processing** – build utilities that consume `SolveResult` to compute derived
  quantities or to export field data to VTK/CSV formats.

## Supporting utilities
Although AlmondFEM is header-only, the repository includes reusable helpers that underpin
the demos and can be leveraged in downstream projects.

### `safe_io`
- Location: `shared/safe_io/include/safe_io/utils.hpp`.
- Purpose: type-safe formatting wrappers around `{fmt}` for consistent logging.
- Integration: link the `safe_io` static library and include `<safe_io/utils.hpp>`.

### Sample applications
- **`cmakeapp1`** – reference CLI application exercising the solver and printing convergence
  diagnostics. Ideal for cross-platform CMake/Ninja workflows.
- **`Application1`** – Visual Studio counterpart that mirrors the CLI behaviour for `.sln`
  users and MSVC debugging sessions.

## Maintenance checklist
- Update this overview whenever new solver headers, utility libraries, or extension hooks
  are added.
- Cross-link new documentation or tutorials from here so contributors can find examples
  quickly.
