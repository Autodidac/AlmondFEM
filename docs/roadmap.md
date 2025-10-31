# AlmondFEM roadmap

The roadmap groups upcoming solver work, tooling refinements, and documentation initiatives
by horizon. Use it to coordinate contributions and identify high-impact areas.

## Short-term (next 1–2 releases)
- **Solver validation suite** – add doctest/catch2-based smoke tests that run sample meshes
  through both direct and CG paths via CMake `ctest`.
- **Build tooling polish** – expand the shell scripts with verbose logging flags and publish
  worked examples for integrating with GitHub Actions.
- **Preset curation** – document each `CMakePresets.json` entry directly in the file and align
  naming across scripts and IDE configurations.
- **Documentation uplift** – refresh API snippets and add a quick-start mesh example to the
  README.
- **Hybrid FLIP/Bubble showcase** – iterate on the sparse particle demo and capture frame
  dumps suitable for inclusion alongside the PDFs.

## Mid-term (next quarter)
- **Preconditioner experiments** – prototype ILU/SSOR variants and surface them through
  `SolverOptions` for comparative studies.
- **Mesh import/export** – ship OBJ/VTK adapters and document integration with the reference
  applications.
- **Performance instrumentation** – provide optional timing hooks and benchmark harnesses that
  can be toggled during configuration.
- **CI coverage** – stand up Linux and Windows pipelines that exercise GCC, Clang, and MSVC
  builds, running the validation suite.

## Long-term (6–12 months)
- **3D element support** – extend `Mesh` and assembly routines to handle tetrahedral meshes
  while preserving backward compatibility with 2D workflows.
- **Domain-specific visualisation** – integrate lightweight viewers or exporters (e.g. VTK,
  ParaView scripts) for post-processing fields.
- **Ecosystem integration** – explore Python bindings and notebook-friendly APIs for rapid
  prototyping.
- **Sustainability processes** – establish release cadences, backward compatibility policies,
  and automated changelog generation.

## Contributing to the roadmap
- Discuss proposed updates in issues or discussions before editing the roadmap.
- When landing new features, update the relevant horizon and cross-link supporting
  documentation.
- Periodically revisit horizons to retire completed items and reprioritise upcoming work.
