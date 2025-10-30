#pragma once

#include "detail/dense_matrix.hpp"
#include "mesh.hpp"
#include "problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <safe_io/utils.hpp>

namespace almond::fem
{
    namespace detail
    {
        inline double element_area(const Mesh& mesh, const Element& element)
        {
            const auto& a = mesh.node(element.node_ids[0]);
            const auto& b = mesh.node(element.node_ids[1]);
            const auto& c = mesh.node(element.node_ids[2]);

            const double area = 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
            return std::abs(area);
        }

        inline std::array<double, 9> local_stiffness(const Mesh& mesh, const Element& element)
        {
            const auto& n0 = mesh.node(element.node_ids[0]);
            const auto& n1 = mesh.node(element.node_ids[1]);
            const auto& n2 = mesh.node(element.node_ids[2]);

            const std::array<double, 3> b{
                n1.y - n2.y,
                n2.y - n0.y,
                n0.y - n1.y,
            };

            const std::array<double, 3> c{
                n2.x - n1.x,
                n0.x - n2.x,
                n1.x - n0.x,
            };

            const double area = element_area(mesh, element);
            const double factor = element.conductivity / (4.0 * area);

            std::array<double, 9> k{};
            for (std::size_t i = 0; i < 3; ++i)
            {
                for (std::size_t j = 0; j < 3; ++j)
                {
                    k[i * 3 + j] = factor * (b[i] * b[j] + c[i] * c[j]);
                }
            }
            return k;
        }
    } // namespace detail

    inline SolveResult solve(const Mesh& mesh, const ProblemDefinition& problem, const SolveOptions& options = {})
    {
        mesh.validate();

        const std::size_t node_count = mesh.node_count();
        detail::DenseMatrix global(node_count);
        std::vector<double> rhs(node_count, 0.0);

        for (const auto& element : mesh.elements())
        {
            const auto local = detail::local_stiffness(mesh, element);
            const double area = detail::element_area(mesh, element);

            for (std::size_t local_i = 0; local_i < 3; ++local_i)
            {
                const auto global_i = element.node_ids[local_i];

                rhs[global_i] += problem.uniform_source * area / 3.0;

                for (std::size_t local_j = 0; local_j < 3; ++local_j)
                {
                    const auto global_j = element.node_ids[local_j];
                    global(global_i, global_j) += local[local_i * 3 + local_j];
                }
            }
        }

        for (const auto& load : problem.point_loads)
        {
            if (load.node >= node_count)
            {
                throw std::out_of_range("Point load references an invalid node");
            }
            rhs[load.node] += load.value;
        }

        for (const auto& condition : problem.dirichlet_conditions)
        {
            if (condition.node >= node_count)
            {
                throw std::out_of_range("Dirichlet boundary references an invalid node");
            }

            for (std::size_t column = 0; column < node_count; ++column)
            {
                global(condition.node, column) = 0.0;
            }
            for (std::size_t row = 0; row < node_count; ++row)
            {
                global(row, condition.node) = 0.0;
            }

            global(condition.node, condition.node) = 1.0;
            rhs[condition.node] = condition.value;
        }

        if (options.verbose)
        {
            safe_io::print("Assembled global stiffness matrix ({}x{})", node_count, node_count);
            for (std::size_t row = 0; row < node_count; ++row)
            {
                std::string row_values;
                row_values.reserve(node_count * 10);
                for (std::size_t column = 0; column < node_count; ++column)
                {
                    row_values += fmt::format("{:>10.4f}", global(row, column));
                }
                safe_io::print("{}", row_values);
            }
            safe_io::print("RHS vector: {}", fmt::join(rhs, ", "));
        }

        auto solution = detail::solve(global, rhs, options.pivot_tolerance);

        std::vector<double> residual(node_count, 0.0);
        for (std::size_t i = 0; i < node_count; ++i)
        {
            double sum = 0.0;
            for (std::size_t j = 0; j < node_count; ++j)
            {
                sum += global(i, j) * solution[j];
            }
            residual[i] = sum - rhs[i];
        }

        const double residual_norm = std::sqrt(std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0));

        if (options.verbose)
        {
            safe_io::print("Residual L2 norm: {:.6e}", residual_norm);
        }

        return SolveResult{std::move(solution), residual_norm};
    }
} // namespace almond::fem

