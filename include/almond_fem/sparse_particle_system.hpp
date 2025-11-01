
#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

namespace almond::fem::particles
{
    struct Particle
    {
        enum class Phase
        {
            Carrier,
            Bubble
        };

        std::array<double, 3> position{0.0, 0.0, 0.0};
        std::array<double, 3> velocity{0.0, 0.0, 0.0};
        double radius{0.02};
        double mass{1.0};
        Phase phase{Phase::Carrier};
        std::uint32_t bubble_id{std::numeric_limits<std::uint32_t>::max()};
    };

    namespace detail
    {
        inline std::array<double, 3> add(const std::array<double, 3>& a, const std::array<double, 3>& b)
        {
            return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
        }

        inline std::array<double, 3> subtract(const std::array<double, 3>& a, const std::array<double, 3>& b)
        {
            return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
        }

        inline std::array<double, 3> scale(const std::array<double, 3>& v, double s)
        {
            return {v[0] * s, v[1] * s, v[2] * s};
        }

        inline double dot(const std::array<double, 3>& a, const std::array<double, 3>& b)
        {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }

        inline double length(const std::array<double, 3>& v)
        {
            return std::sqrt(dot(v, v));
        }

        inline std::array<double, 3> normalize(const std::array<double, 3>& v)
        {
            const auto len = length(v);
            if (len < 1e-9)
            {
                return {0.0, 1.0, 0.0};
            }
            return scale(v, 1.0 / len);
        }

        inline std::array<double, 3> cross(const std::array<double, 3>& a, const std::array<double, 3>& b)
        {
            return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
        }

        inline std::array<double, 3> orthonormal(const std::array<double, 3>& v)
        {
            const auto axis = (std::abs(v[0]) < 0.5) ? std::array<double, 3>{1.0, 0.0, 0.0}
                                                     : std::array<double, 3>{0.0, 1.0, 0.0};
            return normalize(cross(v, axis));
        }

        class AccumulatingScopeTimer
        {
        public:
            explicit AccumulatingScopeTimer(double& accumulator)
                : accumulator_{&accumulator},
                  start_{Clock::now()}
            {
            }

            AccumulatingScopeTimer(const AccumulatingScopeTimer&) = delete;
            AccumulatingScopeTimer& operator=(const AccumulatingScopeTimer&) = delete;

            ~AccumulatingScopeTimer()
            {
                if (accumulator_ != nullptr)
                {
                    const auto end = Clock::now();
                    const std::chrono::duration<double> elapsed = end - start_;
                    *accumulator_ += elapsed.count();
                }
            }

        private:
            using Clock = std::chrono::steady_clock;

            double* accumulator_{nullptr};
            Clock::time_point start_{};
        };

    } // namespace detail

    class SparseParticleSystem
    {
    public:
        struct StepConfig
        {
            double dt{1.0 / 60.0};
            std::array<double, 3> gravity{0.0, -9.81, 0.0};
            double carrier_drag{1.2};
            double bubble_drag{0.35};
            double buoyancy{4.5};
            double max_speed{12.0};
        };

        struct StepStats
        {
            double max_speed{0.0};
            double average_density{0.0};
            std::size_t active_cells{0};
            double velocity_update_seconds{0.0};
            double position_integration_seconds{0.0};
            double spatial_index_rebuild_seconds{0.0};
            double overlap_resolution_seconds{0.0};
            double density_pass_seconds{0.0};
        };

        SparseParticleSystem(double cell_size,
                             std::array<double, 3> domain_min,
                             std::array<double, 3> domain_max)
            : cell_size_{cell_size},
              domain_min_{domain_min},
              domain_max_{domain_max}
        {
        }

        std::size_t add_particle(const Particle& particle)
        {
            particles_.push_back(particle);
            spatial_dirty_ = true;
            return particles_.size() - 1;
        }

        std::vector<Particle>& particles()
        {
            return particles_;
        }

        const std::vector<Particle>& particles() const
        {
            return particles_;
        }

        Particle& particle(std::size_t index)
        {
            return particles_.at(index);
        }

        const Particle& particle(std::size_t index) const
        {
            return particles_.at(index);
        }

        void clear()
        {
            particles_.clear();
            grid_.clear();
            spatial_dirty_ = false;
        }

        StepStats step(const StepConfig& cfg)
        {
            StepStats stats{};

            {
                detail::AccumulatingScopeTimer timer(stats.spatial_index_rebuild_seconds);
                rebuild_spatial_index();
            }

            std::vector<std::vector<std::size_t>> neighbor_cache;
            const bool cache_neighbors = particles_.size() > 1;
            if (cache_neighbors)
            {
                neighbor_cache.resize(particles_.size());
                for (std::size_t index = 0; index < particles_.size(); ++index)
                {
                    neighbor_cache[index] = neighbor_indices(index);
                }
            }

            {
                detail::AccumulatingScopeTimer timer(stats.velocity_update_seconds);
                for (auto& particle : particles_)
                {
                    const auto drag = (particle.phase == Particle::Phase::Carrier) ? cfg.carrier_drag : cfg.bubble_drag;
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        particle.velocity[axis] += cfg.gravity[axis] * cfg.dt;
                    }

                    if (particle.phase == Particle::Phase::Bubble)
                    {
                        particle.velocity[1] += cfg.buoyancy * cfg.dt;
                    }

                    const auto attenuation = std::max(0.0, 1.0 - drag * cfg.dt);
                    for (auto& component : particle.velocity)
                    {
                        component *= attenuation;
                    }

                    const auto speed = detail::length(particle.velocity);
                    if (speed > cfg.max_speed)
                    {
                        const auto capped = detail::scale(detail::normalize(particle.velocity), cfg.max_speed);
                        particle.velocity = capped;
                    }

                    stats.max_speed = std::max(stats.max_speed, detail::length(particle.velocity));
                }
            }

            {
                detail::AccumulatingScopeTimer timer(stats.position_integration_seconds);
                for (auto& particle : particles_)
                {
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        particle.position[axis] += particle.velocity[axis] * cfg.dt;

                        if (particle.position[axis] < domain_min_[axis])
                        {
                            particle.position[axis] = domain_min_[axis];
                            particle.velocity[axis] *= -0.35;
                        }
                        else if (particle.position[axis] > domain_max_[axis])
                        {
                            particle.position[axis] = domain_max_[axis];
                            particle.velocity[axis] *= -0.35;
                        }
                    }
                }
            }

            spatial_dirty_ = !particles_.empty();

            {
                detail::AccumulatingScopeTimer timer(stats.overlap_resolution_seconds);
                const auto* neighbors_ptr = cache_neighbors ? &neighbor_cache : nullptr;
                resolve_overlaps(cfg.dt, neighbors_ptr);
            }

            {
                detail::AccumulatingScopeTimer timer(stats.spatial_index_rebuild_seconds);
                rebuild_spatial_index();
            }

            stats.active_cells = grid_.size();

            const double cell_volume = cell_size_ * cell_size_ * cell_size_;
            {
                detail::AccumulatingScopeTimer timer(stats.density_pass_seconds);
                double accumulated_density = 0.0;
                for (const auto& entry : grid_)
                {
                    accumulated_density += entry.second.total_mass / cell_volume;
                }
                if (!grid_.empty())
                {
                    stats.average_density = accumulated_density / static_cast<double>(grid_.size());
                }
            }

            return stats;
        }

        std::size_t count(Particle::Phase phase) const
        {
            return static_cast<std::size_t>(std::count_if(particles_.begin(), particles_.end(), [phase](const Particle& particle) {
                return particle.phase == phase;
            }));
        }

        std::array<double, 3> domain_min() const
        {
            return domain_min_;
        }

        std::array<double, 3> domain_max() const
        {
            return domain_max_;
        }

    private:
        struct CellKey
        {
            int x{0};
            int y{0};
            int z{0};

            bool operator==(const CellKey& other) const noexcept
            {
                return x == other.x && y == other.y && z == other.z;
            }
        };

        struct CellKeyHasher
        {
            std::size_t operator()(const CellKey& key) const noexcept
            {
                const auto mix = static_cast<std::size_t>(73856093 * key.x ^ 19349663 * key.y ^ 83492791 * key.z);
                return mix;
            }
        };

        struct CellBucket
        {
            std::vector<std::size_t> indices;
            double total_mass{0.0};
        };

        void rebuild_spatial_index()
        {
            if (!spatial_dirty_)
            {
                return;
            }

            grid_.clear();
            grid_.reserve(particles_.size());

            for (std::size_t index = 0; index < particles_.size(); ++index)
            {
                insert_into_cell(index);
            }

            spatial_dirty_ = false;
        }

        void insert_into_cell(std::size_t index)
        {
            const auto key = cell_for_position(particles_[index].position);
            auto& bucket = grid_[key];
            bucket.indices.push_back(index);
            bucket.total_mass += particles_[index].mass;
        }

        CellKey cell_for_position(const std::array<double, 3>& position) const
        {
            CellKey key{};
            const auto clamp_axis = [&](int axis) {
                const auto clamped = std::clamp(position[axis], domain_min_[axis], domain_max_[axis]);
                const auto relative = (clamped - domain_min_[axis]) / cell_size_;
                return static_cast<int>(std::floor(relative));
            };

            key.x = clamp_axis(0);
            key.y = clamp_axis(1);
            key.z = clamp_axis(2);
            return key;
        }

        std::vector<std::size_t> neighbor_indices(std::size_t index) const
        {
            const auto key = cell_for_position(particles_[index].position);
            std::vector<std::size_t> neighbors;
            neighbors.reserve(32);

            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        CellKey neighbor_key{key.x + dx, key.y + dy, key.z + dz};
                        const auto bucket_it = grid_.find(neighbor_key);
                        if (bucket_it == grid_.end())
                        {
                            continue;
                        }
                        neighbors.insert(neighbors.end(), bucket_it->second.indices.begin(), bucket_it->second.indices.end());
                    }
                }
            }

            return neighbors;
        }

        void resolve_overlaps(double dt,
                              const std::vector<std::vector<std::size_t>>* cached_neighbors = nullptr)
        {
            const double inverse_dt = (dt > 1e-6) ? 1.0 / dt : 0.0;
            const bool use_cached_neighbors = cached_neighbors != nullptr &&
                                               cached_neighbors->size() == particles_.size();
            bool any_correction = false;

            for (std::size_t index = 0; index < particles_.size(); ++index)
            {
                auto neighbors = use_cached_neighbors ? (*cached_neighbors)[index] : neighbor_indices(index);
                auto& current = particles_[index];

                for (auto neighbor_index : neighbors)
                {
                    if (neighbor_index <= index)
                    {
                        continue;
                    }

                    auto& neighbor = particles_[neighbor_index];
                    const auto delta = detail::subtract(neighbor.position, current.position);
                    const auto distance = detail::length(delta);
                    const auto min_distance = current.radius + neighbor.radius;

                    if (distance < min_distance && distance > 1e-9)
                    {
                        const auto penetration = min_distance - distance;
                        const auto direction = detail::scale(delta, 1.0 / distance);
                        const auto correction = detail::scale(direction, 0.5 * penetration);

                        for (int axis = 0; axis < 3; ++axis)
                        {
                            current.position[axis] -= correction[axis];
                            neighbor.position[axis] += correction[axis];

                            current.velocity[axis] -= correction[axis] * inverse_dt * 0.5;
                            neighbor.velocity[axis] += correction[axis] * inverse_dt * 0.5;
                        }

                        any_correction = true;
                    }
                }
            }

            if (any_correction)
            {
                spatial_dirty_ = true;
            }
        }

        double cell_size_{0.1};
        std::array<double, 3> domain_min_{0.0, 0.0, 0.0};
        std::array<double, 3> domain_max_{1.0, 1.0, 1.0};

        std::vector<Particle> particles_{};
        std::unordered_map<CellKey, CellBucket, CellKeyHasher> grid_{};
        bool spatial_dirty_{true};
    };

    class BubbleSolver
    {
    public:
        struct CouplingParameters
        {
            double shell_stiffness{6.0};
            double center_smoothing{1.75};
            double buoyancy{5.0};
        };

        struct Telemetry
        {
            std::size_t active_bubbles{0};
            double mean_radius{0.0};
            double mean_height{0.0};
        };

        BubbleSolver()
            : BubbleSolver(CouplingParameters{})
        {
        }

        explicit BubbleSolver(CouplingParameters parameters)
            : parameters_{parameters}
        {
        }

        std::uint32_t create_bubble(double target_radius)
        {
            Bubble bubble{};
            bubble.id = next_identifier_++;
            bubble.target_radius = target_radius;
            bubbles_.push_back(std::move(bubble));
            return bubbles_.back().id;
        }

        void attach_particles(std::uint32_t bubble_id,
                              const std::vector<std::size_t>& indices,
                              SparseParticleSystem& system)
        {
            auto* bubble = find_bubble(bubble_id);
            if (bubble == nullptr)
            {
                return;
            }

            bubble->particle_indices = indices;

            for (auto index : indices)
            {
                auto& particle = system.particle(index);
                particle.phase = Particle::Phase::Bubble;
                particle.bubble_id = bubble_id;
            }
        }

        void apply_constraints(SparseParticleSystem& system, double dt)
        {
            for (auto& bubble : bubbles_)
            {
                if (bubble.particle_indices.empty())
                {
                    continue;
                }

                std::array<double, 3> center{0.0, 0.0, 0.0};
                double count = 0.0;

                for (auto index : bubble.particle_indices)
                {
                    if (index >= system.particles().size())
                    {
                        continue;
                    }

                    const auto& particle = system.particle(index);
                    center = detail::add(center, particle.position);
                    count += 1.0;
                }

                if (count == 0.0)
                {
                    continue;
                }

                center = detail::scale(center, 1.0 / count);

                double mean_radius = 0.0;

                for (auto index : bubble.particle_indices)
                {
                    if (index >= system.particles().size())
                    {
                        continue;
                    }

                    const auto& particle = system.particle(index);
                    const auto offset = detail::subtract(particle.position, center);
                    mean_radius += detail::length(offset);
                }

                mean_radius /= count;
                const auto error = bubble.target_radius - mean_radius;
                const auto radial_correction = parameters_.shell_stiffness * error;

                for (auto index : bubble.particle_indices)
                {
                    if (index >= system.particles().size())
                    {
                        continue;
                    }

                    auto& particle = system.particle(index);
                    auto offset = detail::subtract(particle.position, center);
                    if (detail::length(offset) < 1e-8)
                    {
                        offset = {0.0, 1.0, 0.0};
                    }
                    const auto direction = detail::normalize(offset);

                    for (int axis = 0; axis < 3; ++axis)
                    {
                        particle.velocity[axis] += direction[axis] * radial_correction * dt;
                        particle.velocity[axis] -= offset[axis] * parameters_.center_smoothing * dt * 0.5;
                    }
                }
            }
        }

        void prune_empty(const SparseParticleSystem& system)
        {
            auto matches = [&system](const Bubble& bubble) {
                return std::any_of(bubble.particle_indices.begin(), bubble.particle_indices.end(), [&system, bubble_id = bubble.id](std::size_t index) {
                    if (index >= system.particles().size())
                    {
                        return false;
                    }
                    return system.particle(index).bubble_id == bubble_id;
                });
            };

            bubbles_.erase(std::remove_if(bubbles_.begin(), bubbles_.end(), [&](const Bubble& bubble) {
                                return !matches(bubble);
                            }),
                           bubbles_.end());
        }

        Telemetry gather_telemetry(const SparseParticleSystem& system) const
        {
            Telemetry telemetry{};

            for (const auto& bubble : bubbles_)
            {
                if (bubble.particle_indices.empty())
                {
                    continue;
                }

                std::array<double, 3> center{0.0, 0.0, 0.0};
                double count = 0.0;

                for (auto index : bubble.particle_indices)
                {
                    if (index >= system.particles().size())
                    {
                        continue;
                    }

                    const auto& particle = system.particle(index);
                    center = detail::add(center, particle.position);
                    count += 1.0;
                }

                if (count == 0.0)
                {
                    continue;
                }

                center = detail::scale(center, 1.0 / count);

                double mean_radius = 0.0;
                for (auto index : bubble.particle_indices)
                {
                    if (index >= system.particles().size())
                    {
                        continue;
                    }
                    const auto offset = detail::subtract(system.particle(index).position, center);
                    mean_radius += detail::length(offset);
                }

                mean_radius /= count;
                telemetry.active_bubbles += 1;
                telemetry.mean_radius += mean_radius;
                telemetry.mean_height += center[1];
            }

            if (telemetry.active_bubbles > 0)
            {
                const auto inv = 1.0 / static_cast<double>(telemetry.active_bubbles);
                telemetry.mean_radius *= inv;
                telemetry.mean_height *= inv;
            }

            return telemetry;
        }

        std::size_t bubble_count() const
        {
            return bubbles_.size();
        }

        double buoyancy() const
        {
            return parameters_.buoyancy;
        }

    private:
        struct Bubble
        {
            std::uint32_t id{0};
            double target_radius{0.1};
            std::vector<std::size_t> particle_indices{};
        };

        Bubble* find_bubble(std::uint32_t bubble_id)
        {
            const auto it = std::find_if(bubbles_.begin(), bubbles_.end(), [bubble_id](const Bubble& bubble) {
                return bubble.id == bubble_id;
            });
            return (it == bubbles_.end()) ? nullptr : &(*it);
        }

        const Bubble* find_bubble(std::uint32_t bubble_id) const
        {
            const auto it = std::find_if(bubbles_.begin(), bubbles_.end(), [bubble_id](const Bubble& bubble) {
                return bubble.id == bubble_id;
            });
            return (it == bubbles_.end()) ? nullptr : &(*it);
        }

        CouplingParameters parameters_{};
        std::uint32_t next_identifier_{1};
        std::vector<Bubble> bubbles_{};
    };

    class BubbleEmitter
    {
    public:
        struct Config
        {
            std::array<double, 3> origin{0.5, 0.15, 0.5};
            std::array<double, 3> direction{0.0, 1.0, 0.0};
            double spawn_interval{0.35};
            double min_radius{0.06};
            double max_radius{0.1};
            std::size_t shell_particles{24};
            double launch_speed{0.8};
            double lateral_jitter{0.05};
            double bubble_mass{0.25};
            double stream_spacing{0.08};
        };

        BubbleEmitter()
            : BubbleEmitter(Config{})
        {
        }

        explicit BubbleEmitter(Config config)
            : config_{config},
              radius_distribution_(config.min_radius, config.max_radius),
              jitter_distribution_(-config.lateral_jitter, config.lateral_jitter),
              rng_(5489u)
        {
        }

        bool update(double dt,
                    double /*time*/,
                    SparseParticleSystem& system,
                    BubbleSolver& solver)
        {
            accumulator_ += dt;
            bool emitted = false;

            while (accumulator_ >= config_.spawn_interval)
            {
                accumulator_ -= config_.spawn_interval;

                const auto radius = (config_.min_radius == config_.max_radius)
                                        ? config_.min_radius
                                        : radius_distribution_(rng_);

                auto center = config_.origin;
                const auto direction = detail::normalize(config_.direction);
                center = detail::add(center, detail::scale(direction, config_.stream_spacing * static_cast<double>(emission_count_)));
                center[0] += jitter_distribution_(rng_);
                center[2] += jitter_distribution_(rng_);

                const auto bubble_id = solver.create_bubble(radius);
                std::vector<std::size_t> indices;
                indices.reserve(config_.shell_particles);

                const auto tangent = detail::orthonormal(direction);
                const auto bitangent = detail::normalize(detail::cross(direction, tangent));

                for (std::size_t i = 0; i < config_.shell_particles; ++i)
                {
                    const double phi = (static_cast<double>(i) + 0.5) / static_cast<double>(config_.shell_particles);
                    const double latitude = std::acos(1.0 - 2.0 * phi);
                    const double longitude = std::fmod(static_cast<double>(i) * 2.39996322972865332, 2.0 * 3.14159265358979323846);

                    std::array<double, 3> local{
                        std::cos(longitude) * std::sin(latitude),
                        std::cos(latitude),
                        std::sin(longitude) * std::sin(latitude)};

                    const auto offset = detail::add(
                        detail::add(detail::scale(tangent, local[0]), detail::scale(direction, local[1])),
                        detail::scale(bitangent, local[2]));

                    Particle particle{};
                    particle.position = detail::add(center, detail::scale(offset, radius));
                    particle.radius = radius / 6.0;
                    particle.mass = config_.bubble_mass / static_cast<double>(config_.shell_particles);
                    particle.phase = Particle::Phase::Bubble;
                    particle.velocity = detail::scale(direction, config_.launch_speed);
                    particle.velocity = detail::add(particle.velocity, detail::scale(tangent, jitter_distribution_(rng_) * 0.25));
                    particle.velocity = detail::add(particle.velocity, detail::scale(bitangent, jitter_distribution_(rng_) * 0.25));
                    particle.bubble_id = bubble_id;

                    const auto index = system.add_particle(particle);
                    indices.push_back(index);
                }

                solver.attach_particles(bubble_id, indices, system);
                ++emission_count_;
                emitted = true;
            }

            return emitted;
        }

    private:
        Config config_{};
        double accumulator_{0.0};
        std::size_t emission_count_{0};

        std::uniform_real_distribution<double> radius_distribution_;
        std::uniform_real_distribution<double> jitter_distribution_;
        std::mt19937 rng_;
    };

} // namespace almond::fem::particles
