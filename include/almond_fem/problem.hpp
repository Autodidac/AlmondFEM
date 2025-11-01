#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace almond::fem
{
    namespace detail
    {
        class CsrMatrix;
        class SellCSigmaMatrix;
    }

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
        std::shared_ptr<const detail::CsrMatrix> csr_matrix{};
        std::shared_ptr<const detail::SellCSigmaMatrix> sellc_sigma_matrix{};
    };
} // namespace almond::fem

