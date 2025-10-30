# API overview

The template bundles a small collection of reusable utilities and establishes conventions for adding new libraries. This document summarises what ships today, planned abstractions, and extension points for template consumers.

## AlmondFEM

The AlmondFEM headers live in `include/almond_fem/` and provide a complete, header-only toolkit for assembling and solving 2D
Poisson problems.

### Key headers and types
- **`mesh.hpp`** – exposes the geometric model through `Node`, `Element`, and the `Mesh` container. `Mesh` owns the nodes and
  linear triangular elements, validates connectivity on construction, and offers helpers such as `add_node`, `add_element`,
  and `nodes()`/`elements()` spans for iteration.
- **`problem.hpp`** – defines the loading and constraint primitives: `DirichletBoundary`, `PointLoad`, and
  `ProblemDefinition`, plus the `SolveResult` returned by solver routines.
- **`solver.hpp`** – bundles the solver front-end (`SolverType`, `PreconditionerType`, `SolverOptions`) and implementation
  utilities. Use `SolverOptions` to select between direct Gaussian elimination and the conjugate-gradient path, configure
  tolerance/iteration caps, toggle verbose diagnostics, and request optional SELL-C-sigma views.

Refer back to the [README’s AlmondFEM overview](../README.md#almondfem-library-overview) for a narrative walkthrough and sample
usage that complements this API summary.

### Mesh assembly workflow
1. Create a `Mesh` by either constructing it with node/element vectors or by incrementally calling `add_node`/`add_element`.
   Each `Element` stores three node indices and an optional material conductivity (`conductivity` defaults to `1.0`). Mesh
   validation ensures referenced nodes exist, so populate nodes before pushing elements.
2. Optionally author custom importers by wrapping `Mesh` construction behind your own loader functions—`Mesh` exposes spans for
   read/write access, making it straightforward to integrate procedural generators or on-disk formats.

### Boundary conditions and loads
- Populate `ProblemDefinition::dirichlet_conditions` with `DirichletBoundary` entries (node index + prescribed value) to clamp
  nodal values.
- Add concentrated forces via `ProblemDefinition::point_loads`; each `PointLoad` couples a node index with a scalar load.
- Set `ProblemDefinition::uniform_source` when a domain-wide source term is required.

### Solver and preconditioner selection
- Pick a solver through `SolverOptions::solver` (`SolverType::Direct` or `SolverType::ConjugateGradient`). The iterative path
  honours `SolverOptions::tolerance`, `SolverOptions::max_iterations`, and `SolverOptions::verbose`.
- Choose a preconditioner with `SolverOptions::preconditioner` (`None`, `Jacobi`, or `IncompleteCholesky0`). The IC(0)
  preconditioner requires symmetric positive-definite stiffness matrices; failures emit informative exceptions.
- Advanced knobs:
  - `pivot_tolerance` guards partial pivoting during the direct solve.
  - `build_sellc_sigma`/`sell_chunk_size` request SELL-C-sigma slices alongside the core CSR storage for experimentation or SIMD
    tuning.

### Extension points
- **Custom preconditioners**: derive from `detail::Preconditioner` in `solver.hpp` and inject your implementation into the
  solver dispatch for domain-specific experimentation.
- **Mesh loaders and generators**: leverage `Mesh` spans for efficient conversions from external file formats or procedural
  pipelines.
- **Solver integrations**: the high-level structures (`ProblemDefinition`, `SolverOptions`) are intentionally lightweight, so
  downstream tools can add option fields or adapters without recompiling the core headers.

## Bundled utilities

### `safe_io`
- **Purpose**: provides a type-safe facade around standard output, centralising formatting through `{fmt}`.
- **Key headers**: `shared/safe_io/include/safe_io/utils.hpp` exposes `safe_io::out()` and the templated `safe_io::print()` helper.
- **Linkage**: built as a static library (`safe_io`) and linked into sample applications via CMake (`target_link_libraries(<target> PRIVATE safe_io)`).
- **Usage pattern**: include `<safe_io/utils.hpp>` and call `safe_io::print("Message {}", value);` for buffered, newline-terminated output.

### Static library scaffolding
- **Module**: `StaticLib1` demonstrates how to structure a reusable library with public headers in `StaticLib1/include` and implementation in `StaticLib1/src`.
- **Integration**: its headers are consumed by both `Application1` and `cmakeapp1`, showcasing shared component wiring across different build systems.

## Planned abstractions
- **Configuration layer**: add a lightweight configuration parser (JSON/TOML) with environment overrides and per-target defaults.
- **Logging wrapper**: expand `safe_io` into a logger facade with severity levels, structured context, and sinks for console/files.
- **Platform services**: provide abstractions for file I/O, threading, and GPU context initialisation with unified error handling.
- **Testing helpers**: bundle doctest/catch2 integration with ready-made fixtures and CMake testing presets.

## Extension points
- **New libraries**: mirror `StaticLib1` by creating a directory under `shared/` with `include/` and `src/` subfolders, then register the target in the top-level `CMakeLists.txt`.
- **Application templates**: copy `cmakeapp1` for cross-platform CMake apps or `Application1` for Visual Studio-centric workflows; adjust target linkage as needed.
- **Toolchain hooks**: extend scripts under `cmake/` and `build.sh` to recognise new compilers or build configurations, keeping presets in sync.
- **Documentation**: document new APIs in this file and cross-link any specialised guides under `docs/` so contributors can discover them quickly.

## Maintenance guidelines
- When adding functionality, update this overview to describe new headers, targets, or usage examples.
- Flag experimental APIs with callouts and track their stabilisation in the roadmap and changelog.
- Include code snippets showing best practices for interacting with new abstractions, especially when external dependencies are involved.
