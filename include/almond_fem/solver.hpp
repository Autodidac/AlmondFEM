#pragma once

#include "detail/dense_matrix.hpp"
#include "detail/sparse_matrix.hpp"
#include "mesh.hpp"
#include "problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <safe_io/utils.hpp>

namespace almond::fem
{
    enum class SolverType
    {
        Direct,
        ConjugateGradient,
    };

    enum class PreconditionerType
    {
        None,
        Jacobi,
        IncompleteCholesky0,
    };

    struct SolverOptions
    {
        SolverType solver{SolverType::ConjugateGradient};
        PreconditionerType preconditioner{PreconditionerType::Jacobi};
        double tolerance{1e-8};
        std::size_t max_iterations{1000};
        bool verbose{false};
        double pivot_tolerance{1e-12};
        bool build_sellc_sigma{false};
        std::size_t sell_chunk_size{32};
    };

    inline std::string_view to_string(SolverType type) noexcept
    {
        switch (type)
        {
        case SolverType::Direct:
            return "Direct";
        case SolverType::ConjugateGradient:
            return "ConjugateGradient";
        }
        return "Unknown";
    }

    inline std::string_view to_string(PreconditionerType type) noexcept
    {
        switch (type)
        {
        case PreconditionerType::None:
            return "None";
        case PreconditionerType::Jacobi:
            return "Jacobi";
        case PreconditionerType::IncompleteCholesky0:
            return "IC0";
        }
        return "Unknown";
    }

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

        struct LinearSolveSummary
        {
            std::vector<double> solution{};
            std::size_t iterations{0};
            double achieved_residual{0.0};
        };

        class Preconditioner
        {
        public:
            virtual ~Preconditioner() = default;
            virtual void apply(const std::vector<double>& r, std::vector<double>& z) const = 0;
        };

        class IdentityPreconditioner final : public Preconditioner
        {
        public:
            void apply(const std::vector<double>& r, std::vector<double>& z) const override
            {
                z = r;
            }
        };

        class JacobiPreconditioner final : public Preconditioner
        {
        public:
            explicit JacobiPreconditioner(const CsrMatrix& csr)
                : m_inverse_diagonal(diagonal(csr))
            {
                for (auto& value : m_inverse_diagonal)
                {
                    if (std::abs(value) <= std::numeric_limits<double>::epsilon())
                    {
                        throw std::runtime_error("Jacobi preconditioner requires non-zero diagonal entries");
                    }
                    value = 1.0 / value;
                }
            }

            void apply(const std::vector<double>& r, std::vector<double>& z) const override
            {
                const auto n = m_inverse_diagonal.size();
                z.resize(n);
                for (std::size_t i = 0; i < n; ++i)
                {
                    z[i] = r[i] * m_inverse_diagonal[i];
                }
            }

        private:
            std::vector<double> m_inverse_diagonal{};
        };

        class Ic0Preconditioner final : public Preconditioner
        {
        public:
            explicit Ic0Preconditioner(const CsrMatrix& csr)
            {
                build(csr);
            }

            void apply(const std::vector<double>& r, std::vector<double>& z) const override
            {
                const auto n = m_rows.size();
                if (z.size() != n)
                {
                    z.assign(n, 0.0);
                }

                std::vector<double> y(n, 0.0);
                for (std::size_t row = 0; row < n; ++row)
                {
                    double sum = r[row];
                    double diag = m_diag[row];
                    for (const auto& entry : m_rows[row])
                    {
                        const auto column = static_cast<std::size_t>(entry.first);
                        if (column == row)
                        {
                            diag = entry.second;
                        }
                        else
                        {
                            sum -= entry.second * y[column];
                        }
                    }
                    y[row] = sum / diag;
                }

                for (std::ptrdiff_t row = static_cast<std::ptrdiff_t>(n) - 1; row >= 0; --row)
                {
                    double sum = y[static_cast<std::size_t>(row)];
                    double diag = m_diag[static_cast<std::size_t>(row)];
                    for (const auto& entry : m_upper[static_cast<std::size_t>(row)])
                    {
                        sum -= entry.second * z[static_cast<std::size_t>(entry.first)];
                    }
                    z[static_cast<std::size_t>(row)] = sum / diag;
                }
            }

        private:
            void build(const CsrMatrix& csr)
            {
                const auto n = csr.dimension();
                m_rows.assign(n, {});
                m_upper.assign(n, {});
                m_diag.assign(n, 0.0);

                const auto& row_ptr = csr.row_ptr();
                const auto& col_idx = csr.col_idx();
                const auto& values = csr.values();

                for (std::size_t row = 0; row < n; ++row)
                {
                    const auto begin = static_cast<std::size_t>(row_ptr[row]);
                    const auto end = static_cast<std::size_t>(row_ptr[row + 1]);
                    for (std::size_t idx = begin; idx < end; ++idx)
                    {
                        const auto column = static_cast<std::size_t>(col_idx[idx]);
                        if (column <= row)
                        {
                            m_rows[row].emplace_back(static_cast<int>(column), values[idx]);
                        }
                    }

                    std::sort(m_rows[row].begin(), m_rows[row].end(), [](const auto& lhs, const auto& rhs) {
                        return lhs.first < rhs.first;
                    });

                    const auto diag_iterator = std::find_if(m_rows[row].begin(), m_rows[row].end(), [row](const auto& entry) {
                        return static_cast<std::size_t>(entry.first) == row;
                    });

                    if (diag_iterator == m_rows[row].end())
                    {
                        throw std::runtime_error("IC(0) preconditioner requires explicit diagonal entries");
                    }
                }

                std::vector<std::size_t> diag_index(n, 0);
                for (std::size_t row = 0; row < n; ++row)
                {
                    const auto iterator = std::find_if(m_rows[row].begin(), m_rows[row].end(), [row](const auto& entry) {
                        return static_cast<std::size_t>(entry.first) == row;
                    });
                    diag_index[row] = static_cast<std::size_t>(std::distance(m_rows[row].begin(), iterator));
                }

                for (std::size_t row = 0; row < n; ++row)
                {
                    auto& row_entries = m_rows[row];
                    const auto diag_pos = diag_index[row];

                    for (std::size_t idx = 0; idx < row_entries.size(); ++idx)
                    {
                        auto& entry = row_entries[idx];
                        const auto column = static_cast<std::size_t>(entry.first);
                        if (column == row)
                        {
                            continue;
                        }

                        double sum = entry.second;
                        const auto& column_entries = m_rows[column];
                        std::size_t col_diag = diag_index[column];

                        std::size_t i = 0;
                        std::size_t j = 0;
                        while (i < idx && j < col_diag)
                        {
                            const auto left_column = static_cast<std::size_t>(row_entries[i].first);
                            const auto right_column = static_cast<std::size_t>(column_entries[j].first);

                            if (left_column < column && right_column < column)
                            {
                                if (left_column == right_column)
                                {
                                    sum -= row_entries[i].second * column_entries[j].second;
                                    ++i;
                                    ++j;
                                }
                                else if (left_column < right_column)
                                {
                                    ++i;
                                }
                                else
                                {
                                    ++j;
                                }
                            }
                            else
                            {
                                if (left_column >= column)
                                {
                                    ++j;
                                }
                                else
                                {
                                    ++i;
                                }
                            }
                        }

                        entry.second = sum / column_entries[col_diag].second;
                    }

                    double diag_value = row_entries[diag_pos].second;
                    for (std::size_t idx = 0; idx < diag_pos; ++idx)
                    {
                        const double lij = row_entries[idx].second;
                        diag_value -= lij * lij;
                    }

                    if (diag_value <= 0.0)
                    {
                        throw std::runtime_error("IC(0) factorisation failed: matrix is not SPD");
                    }

                    const double sqrt_value = std::sqrt(diag_value);
                    row_entries[diag_pos].second = sqrt_value;
                    m_diag[row] = sqrt_value;
                }

                for (std::size_t row = 0; row < n; ++row)
                {
                    for (const auto& entry : m_rows[row])
                    {
                        const auto column = static_cast<std::size_t>(entry.first);
                        if (column < row)
                        {
                            m_upper[column].emplace_back(static_cast<int>(row), entry.second);
                        }
                    }
                }
            }

            std::vector<std::vector<std::pair<int, double>>> m_rows{};
            std::vector<std::vector<std::pair<int, double>>> m_upper{};
            std::vector<double> m_diag{};
        };

        inline bool is_structurally_symmetric(const CsrMatrix& csr, double tolerance = 1e-9)
        {
            const auto n = csr.dimension();
            const auto& row_ptr = csr.row_ptr();
            const auto& col_idx = csr.col_idx();
            const auto& values = csr.values();

            for (std::size_t row = 0; row < n; ++row)
            {
                const auto begin = static_cast<std::size_t>(row_ptr[row]);
                const auto end = static_cast<std::size_t>(row_ptr[row + 1]);
                for (std::size_t idx = begin; idx < end; ++idx)
                {
                    const auto column = static_cast<std::size_t>(col_idx[idx]);
                    if (column < row)
                    {
                        continue;
                    }

                    const double value = values[idx];
                    const auto mirror_begin = static_cast<std::size_t>(row_ptr[column]);
                    const auto mirror_end = static_cast<std::size_t>(row_ptr[column + 1]);
                    bool found = false;
                    for (std::size_t mirror = mirror_begin; mirror < mirror_end; ++mirror)
                    {
                        if (static_cast<std::size_t>(col_idx[mirror]) == row)
                        {
                            if (std::abs(values[mirror] - value) > tolerance)
                            {
                                return false;
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        inline std::unique_ptr<Preconditioner> make_preconditioner(const CsrMatrix& csr, PreconditionerType type)
        {
            switch (type)
            {
            case PreconditionerType::None:
                return std::make_unique<IdentityPreconditioner>();
            case PreconditionerType::Jacobi:
                return std::make_unique<JacobiPreconditioner>(csr);
            case PreconditionerType::IncompleteCholesky0:
                if (!is_structurally_symmetric(csr))
                {
                    throw std::invalid_argument("IC(0) preconditioner requires a symmetric matrix");
                }
                return std::make_unique<Ic0Preconditioner>(csr);
            }

            throw std::invalid_argument("Unsupported preconditioner type");
        }

        inline LinearSolveSummary direct_solve(const CsrMatrix& csr, const std::vector<double>& rhs, const SolverOptions& options)
        {
            LinearSolveSummary summary{};
            summary.solution = solve(csr.to_dense(), rhs, options.pivot_tolerance);
            summary.iterations = 1;
            return summary;
        }

        inline LinearSolveSummary conjugate_gradient(const CsrMatrix& csr, const std::vector<double>& rhs, const SolverOptions& options, const Preconditioner& preconditioner)
        {
            const auto n = csr.dimension();
            LinearSolveSummary summary{};
            summary.solution.assign(n, 0.0);

            std::vector<double> r = rhs;
            std::vector<double> z(n, 0.0);
            preconditioner.apply(r, z);
            std::vector<double> p = z;
            std::vector<double> Ap(n, 0.0);

            double rho = std::inner_product(r.begin(), r.end(), z.begin(), 0.0);
            if (std::abs(rho) <= std::numeric_limits<double>::epsilon())
            {
                summary.achieved_residual = std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
                return summary;
            }

            const double tolerance = options.tolerance > 0.0 ? options.tolerance : 1e-12;
            const std::size_t max_iterations = options.max_iterations != 0 ? options.max_iterations : n * 10;

            double residual_norm = std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
            if (residual_norm < tolerance)
            {
                summary.achieved_residual = residual_norm;
                return summary;
            }

            for (std::size_t iteration = 0; iteration < max_iterations; ++iteration)
            {
                multiply(csr, p, Ap);
                double denom = std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);
                if (std::abs(denom) <= std::numeric_limits<double>::epsilon())
                {
                    throw std::runtime_error("CG breakdown: encountered zero denominator");
                }

                const double alpha = rho / denom;
                for (std::size_t i = 0; i < n; ++i)
                {
                    summary.solution[i] += alpha * p[i];
                    r[i] -= alpha * Ap[i];
                }

                residual_norm = std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
                summary.iterations = iteration + 1;
                if (residual_norm < tolerance)
                {
                    break;
                }

                preconditioner.apply(r, z);
                const double rho_new = std::inner_product(r.begin(), r.end(), z.begin(), 0.0);
                if (std::abs(rho_new) <= std::numeric_limits<double>::epsilon())
                {
                    throw std::runtime_error("CG breakdown: encountered zero rho");
                }
                const double beta = rho_new / rho;
                rho = rho_new;
                for (std::size_t i = 0; i < n; ++i)
                {
                    p[i] = z[i] + beta * p[i];
                }
            }

            summary.achieved_residual = residual_norm;
            return summary;
        }

        inline LinearSolveSummary solve_linear_system(const CsrMatrix& csr, const std::vector<double>& rhs, const SolverOptions& options)
        {
            LinearSolveSummary summary{};

            switch (options.solver)
            {
            case SolverType::Direct:
                summary = direct_solve(csr, rhs, options);
                break;
            case SolverType::ConjugateGradient:
            {
                auto preconditioner = make_preconditioner(csr, options.preconditioner);
                summary = conjugate_gradient(csr, rhs, options, *preconditioner);
                break;
            }
            default:
                throw std::invalid_argument("Unsupported solver type");
            }

            std::vector<double> residual(rhs.size(), 0.0);
            multiply(csr, summary.solution, residual);
            for (std::size_t i = 0; i < residual.size(); ++i)
            {
                residual[i] -= rhs[i];
            }
            summary.achieved_residual = std::sqrt(std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0));

            return summary;
        }
    } // namespace detail

    inline SolveResult solve(const Mesh& mesh, const ProblemDefinition& problem, const SolverOptions& options = {})
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
            safe_io::print("Invoking {} solver with {} preconditioner", to_string(options.solver), to_string(options.preconditioner));
        }

        auto summary = detail::solve_linear_system(csr, rhs, options);

        if (options.verbose)
        {
            safe_io::print("Solver completed in {} iteration(s). Residual L2 norm: {:.6e}", summary.iterations, summary.achieved_residual);
        }

        return SolveResult{std::move(summary.solution), summary.achieved_residual};
    }
} // namespace almond::fem

