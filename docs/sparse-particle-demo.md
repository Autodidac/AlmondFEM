
# Sparse particle bubble demo

This walkthrough explains how the repository's new sparse particle domain pairs with the
`UniBubbles`-style emission controls referenced in the research PDFs bundled under
`docs/`. The goal is to mirror the workflow described in **Adaptive Phase-Field FLIP**
while keeping everything runnable from a single CLI.

## Overview
- The simulation volume is discretised into coarse voxels (`SparseParticleSystem`).
- Dense liquid parcels are seeded at the bottom of the tank and evolve as "carrier" phase
  particles.
- `BubbleEmitter` mirrors the sparse particle seeding strategy from the PDFs by sampling a
  Fibonacci sphere to construct closed bubble shells.
- `BubbleSolver` maintains cohesion of each shell and adds a controllable buoyancy term to
  approximate the air/liquid coupling discussed in the papers.

## Running the demo
1. Configure and build the project as usual (`cmake --preset ninja-release && cmake --build --preset ninja-release`).
2. Execute `cmakeapp1`. The program will stream frame-by-frame summaries that list cell
   occupancy, density, and aggregate bubble radii.
3. Tweak the parameters inside `cmakeapp1/src/main.cpp` to explore different emission
   rates, coupling strengths, or voxel resolutions.

## Parameter map
| Component | Paper inspiration | Knob in code | Effect |
|-----------|------------------|--------------|--------|
| Sparse voxel size | Adaptive particle budgets (Adaptive Phase-Field FLIP) | `SparseParticleSystem` constructor | Controls how many parcels occupy a single grid cell. |
| Bubble shell stiffness | Phase-field pressure solve | `BubbleSolver::CouplingParameters::shell_stiffness` | Tightens or loosens the bubble membrane. |
| Center smoothing | UniBubbles cohesion term | `BubbleSolver::CouplingParameters::center_smoothing` | Dampens wobble and keeps the bubble centered. |
| Buoyancy | Air/liquid density contrast | `BubbleSolver::CouplingParameters::buoyancy` | Raises or lowers bubbles relative to gravity. |
| Emission interval | Sparse FLIP reseeding cadence | `BubbleEmitter::Config::spawn_interval` | How frequently new bubbles appear. |
| Shell particles | Sparse shell resolution | `BubbleEmitter::Config::shell_particles` | Number of surface samples per bubble. |

## Extending further
- Couple the particle system back into the FEM solver by sampling bubble pressure into the
  mesh and re-solving for structural responses.
- Export the per-step telemetry to a CSV or VTK file for offline visualisation in ParaView.
- Replace the ASCII driver with a small renderer (SDL, ImGui, bgfx) if you need visual
  footage of the coupled simulation.
