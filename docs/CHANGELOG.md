# Changelog

All notable changes to AlmondFEM will be documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Sparse matrix backends (CSR, optional SELL-C-sigma slices, and an ELLPACK reference
  implementation) selectable through `SolverOptions`.

### Changed
- `SolverOptions` exposes toggles for opting into SELL-C-sigma construction while keeping CSR
  the canonical CPU format.

### Fixed
- Document build-system behaviour in the README and supporting guides.

### Removed
- _None yet._ Use this section to track deprecations.

## Versioning strategy
- Releases follow `MAJOR.MINOR.PATCH` semantics.
- Increment **MAJOR** when making incompatible API or tooling changes.
- Increment **MINOR** when adding functionality in a backward compatible manner.
- Increment **PATCH** for backward compatible bug fixes or documentation updates.
- Each release entry should include the release date and link to comparison diffs when
  available.

## Maintenance checklist
- Update the `[Unreleased]` section as changes merge into the main branch.
- When cutting a release, create a new section for the version (e.g., `## [1.2.0] - 2024-05-14`)
  and reset `[Unreleased]` with placeholders.
- Ensure the README and roadmap reflect any notable changes introduced in the release.
