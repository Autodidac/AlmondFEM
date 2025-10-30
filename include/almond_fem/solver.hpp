#pragma once

#include "detail/dense_matrix.hpp"
#include "detail/sparse_matrix.hpp"
#include "mesh.hpp"
#include "problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <unordered_map>

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
        std::vector<double> rhs(node_count, 0.0);

        std::vector<bool> is_dirichlet(node_count, false);
        std::vector<double> dirichlet_values(node_count, 0.0);
        for (const auto& condition : problem.dirichlet_conditions)
        {
            if (condition.node >= node_count)
            {
                throw std::out_of_range("Dirichlet boundary references an invalid node");
            }
            is_dirichlet[condition.node] = true;
            dirichlet_values[condition.node] = condition.value;
        }

        std::vector<std::vector<std::size_t>> adjacency(node_count);
        for (const auto& element : mesh.elements())
        {
            for (std::size_t local_i = 0; local_i < 3; ++local_i)
            {
                const auto global_i = element.node_ids[local_i];
                if (is_dirichlet[global_i])
                {
                    continue;
                }

                auto& row_adjacency = adjacency[global_i];
                for (std::size_t local_j = 0; local_j < 3; ++local_j)
                {
                    const auto global_j = element.node_ids[local_j];
                    if (is_dirichlet[global_j])
                    {
                        continue;
                    }

                    if (std::find(row_adjacency.begin(), row_adjacency.end(), global_j) == row_adjacency.end())
                    {
                        row_adjacency.push_back(global_j);
                    }
                }
            }
        }

        for (std::size_t row = 0; row < node_count; ++row)
        {
            auto& row_adjacency = adjacency[row];
            if (is_dirichlet[row])
            {
                row_adjacency.assign(1, row);
            }
            else
            {
                if (std::find(row_adjacency.begin(), row_adjacency.end(), row) == row_adjacency.end())
                {
                    row_adjacency.push_back(row);
                }
                std::sort(row_adjacency.begin(), row_adjacency.end());
                row_adjacency.erase(std::unique(row_adjacency.begin(), row_adjacency.end()), row_adjacency.end());
            }
        }

        std::vector<int> row_prefix(node_count + 1, 0);
        for (std::size_t row = 0; row < node_count; ++row)
        {
            row_prefix[row + 1] = row_prefix[row] + static_cast<int>(adjacency[row].size());
        }

        detail::CooMatrix coo{};
        coo.dimension = node_count;
        coo.row_prefix = row_prefix;
        const std::size_t total_entries = static_cast<std::size_t>(coo.row_prefix.back());
        coo.rows.assign(total_entries, 0);
        coo.cols.assign(total_entries, 0);
        coo.values.assign(total_entries, 0.0);

        std::vector<std::unordered_map<std::size_t, std::size_t>> index_map(node_count);
        for (std::size_t row = 0; row < node_count; ++row)
        {
            const auto start = static_cast<std::size_t>(coo.row_prefix[row]);
            const auto& columns = adjacency[row];
            index_map[row].reserve(columns.size());
            for (std::size_t offset = 0; offset < columns.size(); ++offset)
            {
                const auto index = start + offset;
                const auto column = columns[offset];
                coo.rows[index] = static_cast<int>(row);
                coo.cols[index] = static_cast<int>(column);
                coo.values[index] = 0.0;
                index_map[row].emplace(column, index);
            }
        }

        for (const auto& element : mesh.elements())
        {
            const auto local = detail::local_stiffness(mesh, element);
            const double area = detail::element_area(mesh, element);

            for (std::size_t local_i = 0; local_i < 3; ++local_i)
            {
                const auto global_i = element.node_ids[local_i];
                if (is_dirichlet[global_i])
                {
                    continue;
                }

                rhs[global_i] += problem.uniform_source * area / 3.0;

                for (std::size_t local_j = 0; local_j < 3; ++local_j)
                {
                    const auto global_j = element.node_ids[local_j];
                    const double contribution = local[local_i * 3 + local_j];

                    if (is_dirichlet[global_j])
                    {
                        rhs[global_i] -= contribution * dirichlet_values[global_j];
                        continue;
                    }

                    const auto iterator = index_map[global_i].find(global_j);
                    if (iterator == index_map[global_i].end())
                    {
                        throw std::logic_error("Missing adjacency entry during assembly");
                    }
                    coo.values[iterator->second] += contribution;
                }
            }
        }

        for (const auto& load : problem.point_loads)
        {
            if (load.node >= node_count)
            {
                throw std::out_of_range("Point load references an invalid node");
            }
            if (!is_dirichlet[load.node])
            {
                rhs[load.node] += load.value;
            }
        }

        for (const auto& condition : problem.dirichlet_conditions)
        {
            rhs[condition.node] = condition.value;
            const auto iterator = index_map[condition.node].find(condition.node);
            if (iterator == index_map[condition.node].end())
            {
                throw std::logic_error("Missing diagonal entry for Dirichlet constraint");
            }
            coo.values[iterator->second] = 1.0;
        }

        detail::CsrMatrix csr(coo);

        if (options.build_sellc_sigma)
        {
            [[maybe_unused]] detail::SellCSigmaMatrix sell(csr, options.sell_chunk_size);
        }

        if (options.verbose)
        {
            safe_io::print("Assembled global stiffness matrix ({}x{})", node_count, node_count);
            const auto dense_view = csr.to_dense();
            for (std::size_t row = 0; row < node_count; ++row)
            {
                std::string row_values;
                row_values.reserve(node_count * 10);
                for (std::size_t column = 0; column < node_count; ++column)
                {
                    row_values += fmt::format("{:>10.4f}", dense_view(row, column));
                }
                safe_io::print("{}", row_values);
            }
            safe_io::print("RHS vector: {}", fmt::join(rhs, ", "));
        }

        auto solution = detail::solve(csr, rhs, options.pivot_tolerance);

        std::vector<double> residual(node_count, 0.0);
        for (std::size_t i = 0; i < node_count; ++i)
        {
            double sum = 0.0;
            const auto row_begin = static_cast<std::size_t>(csr.row_ptr()[i]);
            const auto row_end = static_cast<std::size_t>(csr.row_ptr()[i + 1]);
            for (std::size_t entry = row_begin; entry < row_end; ++entry)
            {
                const auto column = static_cast<std::size_t>(csr.col_idx()[entry]);
                sum += csr.values()[entry] * solution[column];
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

