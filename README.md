# Cpp20_Ultimate_Project_Template

The Ultimate Hello World – cross-platform, multi-editor ready C++20 template for quickly bootstrapping real-world applications and libraries.

## Project vision
- Provide a batteries-included baseline for modern C++20 development that works across Windows and Linux.
- Showcase clean CMake integration with manifest-based dependency management (vcpkg) and graphics-ready defaults (Vulkan SDK).
- Offer repeatable scripts so teams can configure, build, install, and run projects from a single source of truth.
- Ship a modern header-only finite element solver (AlmondFEM) that can be embedded directly into AlmondShell/AlmondEngine tooling.

## AlmondFEM library overview

AlmondFEM lives under `include/almond_fem` and is delivered as a header-only, interface CMake target. It focuses on 2D scalar
Poisson problems with linear triangular elements and is ideal for gameplay prototyping or runtime tooling where a small,
dependency-free solver is desirable.

Key features:

- **C++23 header-only design** – integrate by linking against the `AlmondFEM` interface target; no compilation units are required.
- **Cross-platform** – relies exclusively on the C++ standard library and the bundled `safe_io` logging helpers.
- **Dual solver paths** – select between a direct Gaussian-elimination solve (`SolverType::Direct`) and an iterative conjugate
  gradient path (`SolverType::ConjugateGradient`) from `SolverOptions::solver`.
- **Configurable convergence** – tune `SolverOptions::tolerance`, `SolverOptions::max_iterations`, and
  `SolverOptions::pivot_tolerance` for precise control over accuracy, iteration caps, and pivoting safeguards, with
  `SolverOptions::verbose` exposing per-iteration diagnostics.
- **Preconditioning support** – pair the iterative solver with `SolverOptions::preconditioner` to switch between Jacobi,
  incomplete Cholesky (IC0), or an identity/no-preconditioner configuration.
- **Debug-friendly** – optional verbose mode prints the stiffness matrix, right-hand side, and residual norm via
  `safe_io::print`.

### Sparse matrix backends

AlmondFEM assembles element contributions into a COO structure and promotes it to CSR before dispatching either solver.
CSR remains the authoritative storage throughout factorisation, residual evaluations, and preconditioner builds, with optional
views derived from it when toggled through `SolverOptions`:

- **SELL-C-sigma slices** – enable `SolverOptions::build_sellc_sigma` (customise with `SolverOptions::sell_chunk_size`) to materialise
  SIMD-friendly slices while retaining CSR as the source-of-truth storage.
- **ELL reference layout** – request an ELLPACK view for fixed-width row storage during experimentation; it is generated on
  demand from the CSR assembled under the active `SolverOptions` configuration for validation and documentation scenarios.

```cpp
almond::fem::SolverOptions options{};
options.solver = almond::fem::SolverType::ConjugateGradient; // Showcase the iterative path.
options.preconditioner = almond::fem::PreconditionerType::Jacobi; // Lightweight diagonal scaling.
options.tolerance = 1e-10; // Tighten convergence for the demo mesh.
options.max_iterations = 256; // Guard against runaway solves in documentation builds.
options.build_sellc_sigma = true; // Materialise SELL-C-sigma slices alongside CSR.
options.sell_chunk_size = 16; // Highlight tunable SIMD-friendly chunk widths.
options.verbose = true; // Emit per-iteration residuals to the console.

const auto result = almond::fem::solve(mesh, problem, options);
safe_io::print("Final residual: {:.3e}", result.residual_norm);
```

`Application1` wires up the same configuration so running it assembles the SELL-C-sigma view, enforces the stricter conjugate
gradient tolerance/iteration cap, and prints the verbose residual history alongside the final norm. Deeper solver internals (matrix
assembly, preconditioners, and backend toggles) are documented in [docs/api-overview.md](docs/api-overview.md); combine that guide
with the sample to experiment with new problem definitions or convergence settings. The surrounding template remains easy to extend
to new compilers, targets, and IDE workflows.

## Supported toolchains
| Compiler | Generator | Notes |
| --- | --- | --- |
| GCC | Ninja or Unix Makefiles | Auto-detected by scripts. Works on Linux and Windows (via MSYS2/MinGW). |
| Clang | Ninja or Unix Makefiles | Uses system `clang++`; compatible with LLVM toolchains on all platforms. |
| MSVC | Ninja Multi-Config | Requires Visual Studio Build Tools and an MSVC developer prompt. |

> Configure-time detection automatically injects the appropriate vcpkg toolchain if `VCPKG_ROOT`, `VCPKG_INSTALLATION_ROOT`, or a local `vcpkg_installed` folder is present.

## Quick start
Follow the configure → build → install → run flow:

1. **Configure**
   ```bash
   ./cmake/configure.sh gcc Debug
   ```
   Adjust the compiler (`gcc`, `clang`, `msvc`) and build type (`Debug`, `Release`, etc.) as needed. Configuration output is placed in `Bin/<Compiler>-<Config>`.

2. **Build**
   ```bash
   ./build.sh gcc Debug
   ```
   Reuses the configuration directory if it already exists. Pass `Release` to build optimized binaries.

3. **Install**
   ```bash
   ./install.sh gcc Debug
   ```
   Installs artifacts into `built/bin/<Compiler>-<Config>` so they can be redistributed or consumed by downstream projects.

4. **Run**
   ```bash
   ./run.sh gcc Debug
   ```
   Automatically locates the configured executable (e.g., `cmakeapp1`) and launches it.

Platform-specific setup steps (Visual Studio, VS Code tasks, Linux package managers, etc.) are documented in:
- [docs/windows.md](docs/windows.md)
- [docs/linux.md](docs/linux.md)

Additional high-level documentation lives in:
- [docs/roadmap.md](docs/roadmap.md) for upcoming enhancements and priorities.
- [docs/api-overview.md](docs/api-overview.md) for bundled utilities, abstractions, and extension points.
- [docs/CHANGELOG.md](docs/CHANGELOG.md) for release notes and versioning strategy.

## Environment prerequisites
- **vcpkg**: Required for dependency resolution. Set either `VCPKG_ROOT` or `VCPKG_INSTALLATION_ROOT` before configuring, or place a local bootstrap under `vcpkg_installed/`.
- **Vulkan SDK**: Install the Vulkan SDK for graphics samples or libraries that rely on Vulkan. Ensure `VULKAN_SDK` is exported and its `Bin` directory is in `PATH`.
- **CMake ≥ 3.21**: Needed for multi-config generators and preset support.
- **Ninja (optional)**: Speeds up builds; automatically selected when available.

See the platform guides for package manager recommendations and verification commands.

## Automation scripts
All helper scripts reside at the repository root and are POSIX shell compatible (Windows users should run them inside WSL, MSYS2, or Git Bash).

| Script | Purpose | Usage examples |
| --- | --- | --- |
| `cmake/configure.sh` | Creates the build tree for a compiler + configuration pair. | `./cmake/configure.sh clang Release`<br>`./cmake/configure.sh msvc RelWithDebInfo` |
| `build.sh` | Builds the selected target set. | `./build.sh gcc Debug` |
| `install.sh` | Installs CMake targets into `built/bin/<Compiler>-<Config>`. | `./install.sh clang Release` |
| `run.sh` | Locates and executes the sample executable for the requested configuration. | `./run.sh msvc Release` |

### Troubleshooting
- **Missing vcpkg toolchain**: Ensure `VCPKG_ROOT` or `VCPKG_INSTALLATION_ROOT` points to a valid vcpkg clone. Re-run `./cmake/configure.sh ...` after setting the variable.
- **Generator errors**: Install Ninja (`sudo apt install ninja-build` or `choco install ninja`) or ensure `make` is available. The scripts will fallback to Unix Makefiles when Ninja is absent.
- **MSVC builds fail outside developer prompt**: Use the “x64 Native Tools Command Prompt for VS” so that `cl.exe` and Ninja are on `PATH`.
- **Vulkan loader not found**: Confirm the Vulkan SDK is installed and that `VK_LAYER_PATH` and `VULKAN_SDK` environment variables are exported before configuring.

## Customizing for new targets
1. **Extend the configuration script**: Update `cmake/configure.sh` with a new compiler case (e.g., `icc`, `emcc`) and provide generator/toolchain arguments.
2. **Add CMake presets**: Edit `CMakePresets.json` to encode the new build matrix for IDEs that support presets.
3. **Update automation scripts**: Mirror the new compiler identifier in `build.sh`, `install.sh`, and `run.sh` so the end-to-end flow is consistent.
4. **Document the workflow**: Add a subsection to the appropriate platform guide (`docs/*.md`) describing unique requirements or environment variables.

## Contribution guidelines
- **Roadmap**: When proposing or completing notable work, update [docs/roadmap.md](docs/roadmap.md) to reflect new priorities or delivered items.
- **API overview**: Document new libraries, utilities, or extension points in [docs/api-overview.md](docs/api-overview.md) with usage notes.
- **Changelog**: Record user-visible changes under the `[Unreleased]` section of [docs/CHANGELOG.md](docs/CHANGELOG.md) and promote entries during releases.
- **Cross-links**: Ensure new documentation is linked from the README or relevant guides so contributors can discover it easily.

## Licensing
This template ships with the default license for demonstration purposes. Replace it with your organization’s preferred license before distributing derived work.
