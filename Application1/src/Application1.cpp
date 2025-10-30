#include <almond_fem/almond_fem.hpp>
#include <safe_io/utils.hpp>

#include <array>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace
{
    almond::fem::Mesh build_demo_mesh()
    {
        using almond::fem::Element;
        using almond::fem::Mesh;
        using almond::fem::Node;

        Mesh mesh;
        mesh.add_node(Node{0.0, 0.0});
        mesh.add_node(Node{1.0, 0.0});
        mesh.add_node(Node{1.0, 1.0});
        mesh.add_node(Node{0.0, 1.0});
        mesh.add_node(Node{0.5, 0.5});

        mesh.add_element(Element{{0, 1, 4}, 1.0});
        mesh.add_element(Element{{1, 2, 4}, 1.0});
        mesh.add_element(Element{{2, 3, 4}, 1.0});
        mesh.add_element(Element{{3, 0, 4}, 1.0});

        return mesh;
    }
}

int main()
{
    const auto mesh = build_demo_mesh();

    almond::fem::ProblemDefinition problem{};
    problem.uniform_source = 0.0;
    problem.dirichlet_conditions = {
        {0, 1.0}, // Left boundary: warm
        {3, 1.0},
        {1, 0.0}, // Right boundary: cold
        {2, 0.0},
    };

    // Excite the interior node with a small point load to break symmetry.
    problem.point_loads = {
        {4, 0.1},
    };

    const auto report_solution = [&](std::string_view label, const almond::fem::SolveResult& result) {
        safe_io::print("{} residual norm: {:.6e}", label, result.residual_norm);
        for (std::size_t i = 0; i < result.nodal_values.size(); ++i)
        {
            const auto& node = mesh.node(i);
            safe_io::print("  Node {} @ ({:.2f}, {:.2f}) -> {:.4f}", i, node.x, node.y, result.nodal_values[i]);
        }
    };

    almond::fem::SolverOptions cg_jacobi{};
    cg_jacobi.verbose = true;
    const auto cg_jacobi_result = almond::fem::solve(mesh, problem, cg_jacobi);
    report_solution("CG + Jacobi", cg_jacobi_result);

    almond::fem::SolverOptions cg_ic0{};
    cg_ic0.preconditioner = almond::fem::PreconditionerType::IncompleteCholesky0;
    const auto cg_ic0_result = almond::fem::solve(mesh, problem, cg_ic0);
    report_solution("CG + IC(0)", cg_ic0_result);

    almond::fem::SolverOptions direct{};
    direct.solver = almond::fem::SolverType::Direct;
    direct.preconditioner = almond::fem::PreconditionerType::None;
    const auto direct_result = almond::fem::solve(mesh, problem, direct);
    report_solution("Direct", direct_result);

#ifdef _WIN32
    safe_io::print("Press Enter to exit...");
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
#endif

    return 0;
}
