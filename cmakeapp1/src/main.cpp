
// Demonstration that blends a sparse particle system with buoyant bubbles
#include <almond_fem/sparse_particle_system.hpp>
#include <safe_io/utils.hpp>

#include <chrono>
#include <thread>

namespace
{
    using almond::fem::particles::BubbleEmitter;
    using almond::fem::particles::BubbleSolver;
    using almond::fem::particles::Particle;
    using almond::fem::particles::SparseParticleSystem;

    SparseParticleSystem build_fluid_domain()
    {
        constexpr double cell_size = 0.08;
        SparseParticleSystem system{cell_size, {0.0, 0.0, 0.0}, {1.0, 1.2, 1.0}};

        const auto domain_min = system.domain_min();
        const auto domain_max = system.domain_max();

        const int count_x = 12;
        const int count_y = 6;
        const int count_z = 12;

        const double spacing_x = (domain_max[0] - domain_min[0]) / static_cast<double>(count_x);
        const double spacing_y = (domain_max[1] - domain_min[1]) / static_cast<double>(count_y * 2);
        const double spacing_z = (domain_max[2] - domain_min[2]) / static_cast<double>(count_z);

        for (int ix = 0; ix < count_x; ++ix)
        {
            for (int iy = 0; iy < count_y; ++iy)
            {
                for (int iz = 0; iz < count_z; ++iz)
                {
                    Particle particle{};
                    particle.phase = Particle::Phase::Carrier;
                    particle.mass = 1.5;
                    particle.radius = 0.025;
                    particle.position = {
                        domain_min[0] + (ix + 0.5) * spacing_x,
                        domain_min[1] + (iy + 0.5) * spacing_y,
                        domain_min[2] + (iz + 0.5) * spacing_z};
                    system.add_particle(particle);
                }
            }
        }

        return system;
    }

    void run_demo()
    {
        auto system = build_fluid_domain();

        BubbleSolver solver{BubbleSolver::CouplingParameters{
            .shell_stiffness = 8.5,
            .center_smoothing = 1.1,
            .buoyancy = 7.5}};

        BubbleEmitter emitter{BubbleEmitter::Config{
            .origin = {0.5, 0.08, 0.5},
            .direction = {0.0, 1.0, 0.0},
            .spawn_interval = 0.35,
            .min_radius = 0.05,
            .max_radius = 0.08,
            .shell_particles = 32,
            .launch_speed = 0.95,
            .lateral_jitter = 0.045,
            .bubble_mass = 0.18,
            .stream_spacing = 0.06}};

        SparseParticleSystem::StepConfig step_config{};
        step_config.dt = 1.0 / 90.0;
        step_config.gravity = {0.0, -9.81, 0.0};
        step_config.carrier_drag = 1.6;
        step_config.bubble_drag = 0.28;
        step_config.buoyancy = solver.buoyancy();
        step_config.max_speed = 9.0;

        const std::size_t total_steps = 220;
        double time = 0.0;

        for (std::size_t step = 0; step < total_steps; ++step)
        {
            emitter.update(step_config.dt, time, system, solver);
            solver.apply_constraints(system, step_config.dt);

            const auto stats = system.step(step_config);
            solver.prune_empty(system);
            const auto telemetry = solver.gather_telemetry(system);

            time += step_config.dt;

            if (step % 15 == 0)
            {
                const auto carrier_count = system.count(Particle::Phase::Carrier);
                const auto bubble_count = system.count(Particle::Phase::Bubble);

                safe_io::print(
                    "step {:>3}  t={:.2f}s  carrier={} bubble={} active_cells={} avg_density={:.3f} max|v|={:.3f}",
                    step,
                    time,
                    carrier_count,
                    bubble_count,
                    stats.active_cells,
                    stats.average_density,
                    stats.max_speed);

                if (telemetry.active_bubbles > 0)
                {
                    safe_io::print("          bubble μ radius={:.3f} height={:.3f}", telemetry.mean_radius, telemetry.mean_height);
                }
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }

        safe_io::print("Demo completed after {:.2f}s of simulated time.", time);
    }
} // namespace

int main()
{
    run_demo();
    return 0;
}
