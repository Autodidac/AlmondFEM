# AlmondFEM

AlmondFEM is a modern C++23 finite-element toolkit focused on solving 2D scalar Poisson
problems with linear triangular elements. The repository packages the production headers,
reference applications, and every build system configuration required to exercise the
solver on Linux and Windows. Although the project started from a generic template, the
layout has been curated to serve a simulation-centric workflow: headers live under
`include/`, solver demos sit in `cmakeapp1/` and `Application1/`, and the surrounding
infrastructure targets reproducible scientific builds.

## Project goals
- Deliver a compact, dependency-light solver that can be embedded inside tools or games.
- Preserve cross-platform workflows so researchers can share meshes and executables without
  re-authoring build files.
- Document each build system variant (scripts, presets, IDE integrations, Visual Studio
  solution) so new contributors can pick the flow that matches their environment.
- Provide an extendable foundation for meshing, preconditioning, and linear-algebra
  experiments.

## Solver capabilities
- **Header-only core** – link the `AlmondFEM` interface target to consume the solver without
  adding new translation units.
- **Flexible solver pipeline** – switch between direct Gaussian elimination and conjugate
  gradient iterations via `SolverOptions::solver`.
- **Preconditioning options** – enable Jacobi or incomplete Cholesky (IC(0)) to stabilise
  iterative solves.
- **Optional sparse backends** – materialise SELL-C-sigma or ELLPACK views next to the CSR
  baseline by toggling `SolverOptions`.
- **Diagnostic hooks** – activate `SolverOptions::verbose` to stream assembly output,
  iteration counts, and residual norms through `safe_io`.

Refer to [`docs/api-overview.md`](docs/api-overview.md) for header-by-header details and
usage snippets.

## Repository layout and build systems
AlmondFEM ships multiple build entry points so each platform can use familiar tooling while
sharing the same CMake project graph. The most common paths are summarised below; the
remaining documentation drills into each option.

| Location | What it contains | Typical usage |
| --- | --- | --- |
| `include/almond_fem/` | All solver headers (mesh definitions, problem description, solver front-end). | Consume from your own CMake target or modify the core implementation. |
| `cmakeapp1/` | Cross-platform CMake application that drives AlmondFEM with sample meshes. | Preferred executable for Linux/macOS/WSL; integrates cleanly with presets and Ninja/Make. |
| `Application1/` | Visual Studio oriented application mirroring `cmakeapp1` behaviour. | Native MSVC workflow where `.sln`-based debugging is desired. |
| `shared/` | Supporting libraries such as `safe_io` and reusable utilities. | Shared logging helpers and future runtime extensions. |
| `cmake/` | Toolchain-aware shell scripts and helper modules. | Automates configure/build/install/run steps in a reproducible way. |
| `CMakePresets.json` | IDE and CLI presets for mainstream compilers and configurations. | Allows `cmake --preset` and VS Code/CLion integrations. |
| `app1.sln` | Hand-authored Visual Studio solution. | Directly open in Visual Studio when bypassing presets. |

### Scripted build flow
The top-level shell scripts are the recommended entry point on Linux, macOS, and Windows
(with Git Bash/WSL). They wrap CMake invocations with a consistent directory structure.

1. **Configure** – `./cmake/configure.sh <compiler> <config>` detects the generator, wires
   in vcpkg if available, and initialises `Bin/<Compiler>-<Config>`.
2. **Build** – `./build.sh <compiler> <config>` compiles all default targets; subsequent
   invocations reuse the same build tree.
3. **Install** – `./install.sh <compiler> <config>` stages binaries, libraries, and headers
   under `built/bin/<Compiler>-<Config>` for packaging or downstream consumption.
4. **Run** – `./run.sh <compiler> <config>` locates the demo executable (e.g. `cmakeapp1`)
   and launches it with the chosen toolchain.

Each script prints the exact CMake command it executes, making it easy to adapt for custom
CI pipelines or bespoke compiler flags.

### IDE presets and CMake CLI
`CMakePresets.json` mirrors the script matrix and can be used directly:

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
cmake --build --preset gcc-debug --target install
```

Presets are consumed automatically by Visual Studio, VS Code, and CLion. Add new presets
when introducing compilers, sanitiser builds, or cache wrappers; mirror any additions in
the shell scripts so both pathways stay aligned.

### Visual Studio solution
For developers who rely on `.sln` files, `app1.sln` loads the AlmondFEM library, the
`Application1` executable, and shared utilities. The solution references the same source
files and headers as the CMake targets, ensuring IDE-only work remains compatible with the
scripted and preset flows.

## Quick start
Follow the configure → build → install → run sequence using the compiler for your platform:

```bash
./cmake/configure.sh gcc Debug
./build.sh gcc Debug
./install.sh gcc Debug
./run.sh gcc Debug
```

Swap `gcc` for `clang` or `msvc`, and choose `Release`/`RelWithDebInfo` when benchmarking.
See [`docs/linux.md`](docs/linux.md) and [`docs/windows.md`](docs/windows.md) for platform-
specific package recommendations and IDE guidance.

## Environment prerequisites
- **vcpkg** – optional but recommended for dependency resolution. Export `VCPKG_ROOT` or
  `VCPKG_INSTALLATION_ROOT` before configuring, or drop a local `vcpkg_installed/` folder in
  the repository root.
- **CMake ≥ 3.21** – required for preset support and multi-config generators.
- **Compiler toolchains** – GCC, Clang, or MSVC as outlined in the supported toolchain
  table below.
- **Ninja (optional)** – used automatically when available for faster builds.
- **Vulkan SDK (optional)** – only needed when experimenting with graphics-adjacent demos.

## Supported toolchains
| Compiler | Generator | Notes |
| --- | --- | --- |
| GCC | Ninja or Unix Makefiles | Auto-detected by scripts. Works on Linux and Windows (MSYS2/WSL). |
| Clang | Ninja or Unix Makefiles | Uses system `clang++`; compatible with LLVM toolchains on all platforms. |
| MSVC | Ninja Multi-Config | Requires Visual Studio Build Tools and an MSVC developer prompt. |

> During configuration the scripts attempt to locate a vcpkg toolchain using
> `VCPKG_ROOT`, `VCPKG_INSTALLATION_ROOT`, or a local `vcpkg_installed` directory.

## Additional documentation
- [`docs/api-overview.md`](docs/api-overview.md) – data structures, solver options, and
  extension hooks.
- [`docs/linux.md`](docs/linux.md) – distribution packages, IDE tips, and troubleshooting for
  Unix-like systems.
- [`docs/windows.md`](docs/windows.md) – Windows-specific setup, Visual Studio workflows, and
  environment configuration.
- [`docs/roadmap.md`](docs/roadmap.md) – upcoming solver features, tooling improvements, and
  research milestones.
- [`docs/CHANGELOG.md`](docs/CHANGELOG.md) – release tracking with a Keep-a-Changelog format.
- [`docs/sparse-particle-demo.md`](docs/sparse-particle-demo.md) – step-by-step hybrid sparse
  particle/bubble walkthrough inspired by the FLIP + UniBubbles papers in `docs/`.

## Contributing
- Update the roadmap and changelog whenever solver features or build workflows evolve.
- Document new APIs or utilities under `docs/` and cross-link them from this README.
- Follow the existing directory conventions (`include/`, `shared/`, `docs/`) when introducing
  new modules or samples.
- Prefer adding configuration flags to the scripts/presets instead of hard-coding toolchain
  choices in CMake files.

## Licensing
The repository currently ships with a demonstration license. Replace `LICENSE` with your
organisation's preferred terms before distributing derivatives.
