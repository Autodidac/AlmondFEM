#include <almond_fem/mesh.hpp>
#include <almond_fem/problem.hpp>
#include <almond_fem/solver.hpp>

#include <safe_io/utils.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace
{
    using almond::fem::DirichletBoundary;
    using almond::fem::Mesh;
    using almond::fem::Node;
    using almond::fem::ProblemDefinition;
    using almond::fem::SolveResult;
    using almond::fem::SolverOptions;
    using almond::fem::SolverType;
    using almond::fem::PreconditionerType;

    Mesh build_unit_square_mesh(std::size_t subdivisions)
    {
        const double step = 1.0 / static_cast<double>(subdivisions);
        const std::size_t grid_points = subdivisions + 1;
        std::vector<Node> nodes;
        nodes.reserve(grid_points * grid_points);
        for (std::size_t j = 0; j < grid_points; ++j)
        {
            for (std::size_t i = 0; i < grid_points; ++i)
            {
                nodes.push_back(Node{static_cast<double>(i) * step, static_cast<double>(j) * step});
            }
        }

        std::vector<almond::fem::Element> elements;
        elements.reserve(subdivisions * subdivisions * 2);
        for (std::size_t j = 0; j < subdivisions; ++j)
        {
            for (std::size_t i = 0; i < subdivisions; ++i)
            {
                const std::size_t n0 = j * grid_points + i;
                const std::size_t n1 = n0 + 1;
                const std::size_t n2 = n0 + grid_points;
                const std::size_t n3 = n2 + 1;

                elements.push_back(almond::fem::Element{std::array<std::size_t, 3>{n0, n1, n3}, 1.0});
                elements.push_back(almond::fem::Element{std::array<std::size_t, 3>{n0, n3, n2}, 1.0});
            }
        }

        return Mesh{std::move(nodes), std::move(elements)};
    }

    ProblemDefinition build_problem(const Mesh& mesh)
    {
        ProblemDefinition problem{};
        problem.uniform_source = 1.0;

        const std::size_t node_count = mesh.node_count();
        problem.dirichlet_conditions.reserve(node_count / 4);
        constexpr double eps = 1e-9;
        for (std::size_t idx = 0; idx < node_count; ++idx)
        {
            const auto& node = mesh.node(idx);
            if (std::abs(node.x - 0.0) < eps)
            {
                problem.dirichlet_conditions.push_back(DirichletBoundary{idx, 0.0});
            }
            else if (std::abs(node.x - 1.0) < eps)
            {
                problem.dirichlet_conditions.push_back(DirichletBoundary{idx, 1.0});
            }
            else if (std::abs(node.y - 0.0) < eps || std::abs(node.y - 1.0) < eps)
            {
                problem.dirichlet_conditions.push_back(DirichletBoundary{idx, 0.0});
            }
        }

        return problem;
    }

    struct TimingSummary
    {
        SolveResult result{};
        double average_ms{0.0};
    };

    TimingSummary run_solver(const Mesh& mesh, const ProblemDefinition& problem, bool use_sell, std::size_t repetitions)
    {
        SolverOptions options{};
        options.solver = SolverType::ConjugateGradient;
        options.preconditioner = PreconditionerType::Jacobi;
        options.tolerance = 1e-10;
        options.max_iterations = 4000;
        options.build_sellc_sigma = use_sell;
        options.sell_chunk_size = 32;

        SolveResult last_result{};
        double total_ms = 0.0;

        for (std::size_t iteration = 0; iteration < repetitions; ++iteration)
        {
            const auto start = std::chrono::steady_clock::now();
            last_result = almond::fem::solve(mesh, problem, options);
            const auto end = std::chrono::steady_clock::now();
            total_ms += std::chrono::duration<double, std::milli>(end - start).count();
        }

        TimingSummary summary{};
        summary.result = std::move(last_result);
        summary.average_ms = repetitions > 0 ? total_ms / static_cast<double>(repetitions) : 0.0;
        return summary;
    }

    double max_difference(const std::vector<double>& lhs, const std::vector<double>& rhs)
    {
        if (lhs.size() != rhs.size())
        {
            return std::numeric_limits<double>::infinity();
        }

        double diff = 0.0;
        for (std::size_t i = 0; i < lhs.size(); ++i)
        {
            diff = std::max(diff, std::abs(lhs[i] - rhs[i]));
        }
        return diff;
    }
}

int main()
{
    constexpr std::size_t subdivisions = 40;
    constexpr std::size_t repetitions = 3;

    const auto mesh = build_unit_square_mesh(subdivisions);
    const auto problem = build_problem(mesh);

    const auto csr_run = run_solver(mesh, problem, false, repetitions);
    const auto sell_run = run_solver(mesh, problem, true, repetitions);

    const double nodal_diff = max_difference(csr_run.result.nodal_values, sell_run.result.nodal_values);
    const double residual_diff = std::abs(csr_run.result.residual_norm - sell_run.result.residual_norm);

    safe_io::print(
        "CG benchmark on {}x{} grid ({} nodes, {} elements)",
        subdivisions,
        subdivisions,
        mesh.node_count(),
        mesh.element_count());
    safe_io::print(
        " CSR average {:.3f} ms | residual {:.3e}",
        csr_run.average_ms,
        csr_run.result.residual_norm);
    safe_io::print(
        " SELL-C-sigma average {:.3f} ms | residual {:.3e} | chunk {}",
        sell_run.average_ms,
        sell_run.result.residual_norm,
        sell_run.result.sellc_sigma_matrix ? sell_run.result.sellc_sigma_matrix->chunk_size() : 0);

    if (!sell_run.result.sellc_sigma_matrix)
    {
        safe_io::print("SELL-C-sigma matrix was not retained in the solve result.");
        return 1;
    }

    if (nodal_diff > 1e-8 || residual_diff > 1e-8)
    {
        safe_io::print(
            "Numerical mismatch detected between CSR and SELL results (nodal diff {:.3e}, residual diff {:.3e}).",
            nodal_diff,
            residual_diff);
        return 1;
    }

    safe_io::print(
        "Results match within tolerance. Speedup factor {:.2f}x",
        (csr_run.average_ms > 0.0) ? (csr_run.average_ms / sell_run.average_ms) : 0.0);

    return 0;
}
