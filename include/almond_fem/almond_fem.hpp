#pragma once

#include "mesh.hpp"
#include "problem.hpp"
#include "solver.hpp"
#include "sparse_particle_system.hpp"
#include "bubble_simulation.hpp"

/*
 * AlmondFEM – A lightweight, header-only finite element solver designed for
 * AlmondShell/AlmondEngine integrations.
 *
 * Features:
 *  - 2D scalar Poisson equation support with linear triangular elements.
 *  - Header-only design with zero external dependencies beyond the standard
 *    library and the shared safe_io printing utilities already bundled with
 *    the project.
 *  - Deterministic assembly and dense solver suitable for small to medium
 *    systems such as gameplay prototyping or tooling integrations.
 *  - Cross-platform C++23 interface relying solely on standard containers and
 *    algorithms.
 */

