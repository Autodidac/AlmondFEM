#pragma once

#include <cstddef>
#include <vector>

namespace almond::fem
{
    struct DirichletBoundary
    {
        std::size_t node{};
        double value{0.0};
    };

    struct PointLoad
    {
        std::size_t node{};
        double value{0.0};
    };

    struct SolveOptions
    {
        bool verbose{false};
        double pivot_tolerance{1e-12};
    };

    struct ProblemDefinition
    {
        std::vector<PointLoad> point_loads{};
        std::vector<DirichletBoundary> dirichlet_conditions{};
        double uniform_source{0.0};
    };

    struct SolveResult
    {
        std::vector<double> nodal_values{};
        double residual_norm{0.0};
    };
} // namespace almond::fem

