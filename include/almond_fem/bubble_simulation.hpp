#pragma once

#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace almond::fem::bubbles
{
    struct Vec3
    {
        double x{0.0};
        double y{0.0};
        double z{0.0};

        Vec3& operator+=(const Vec3& other) noexcept
        {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        Vec3& operator-=(const Vec3& other) noexcept
        {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        Vec3& operator*=(double scalar) noexcept
        {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        Vec3& operator/=(double scalar)
        {
            if (std::abs(scalar) <= std::numeric_limits<double>::epsilon())
            {
                throw std::runtime_error("Attempting to divide Vec3 by zero");
            }
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }
    };

    [[nodiscard]] inline Vec3 operator+(Vec3 lhs, const Vec3& rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }

    [[nodiscard]] inline Vec3 operator-(Vec3 lhs, const Vec3& rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }

    [[nodiscard]] inline Vec3 operator*(Vec3 lhs, double scalar) noexcept
    {
        lhs *= scalar;
        return lhs;
    }

    [[nodiscard]] inline Vec3 operator*(double scalar, Vec3 rhs) noexcept
    {
        rhs *= scalar;
        return rhs;
    }

    [[nodiscard]] inline Vec3 operator/(Vec3 lhs, double scalar)
    {
        lhs /= scalar;
        return lhs;
    }

    [[nodiscard]] inline double dot(const Vec3& a, const Vec3& b) noexcept
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    [[nodiscard]] inline double length_squared(const Vec3& value) noexcept
    {
        return dot(value, value);
    }

    [[nodiscard]] inline double length(const Vec3& value) noexcept
    {
        return std::sqrt(length_squared(value));
    }

    [[nodiscard]] inline Vec3 normalise(const Vec3& value)
    {
        const double len = length(value);
        if (len <= std::numeric_limits<double>::epsilon())
        {
            return Vec3{};
        }
        return value / len;
    }

    struct Bubble
    {
        Vec3 position{};
        Vec3 velocity{};
        double radius{0.5};
        double mass{1.0};

        [[nodiscard]] double volume() const noexcept
        {
            return (4.0 / 3.0) * std::numbers::pi * radius * radius * radius;
        }

        [[nodiscard]] double surface_area() const noexcept
        {
            return 4.0 * std::numbers::pi * radius * radius;
        }
    };

    class EulerianGrid
    {
    public:
        EulerianGrid(std::size_t nx, std::size_t ny, std::size_t nz, double cell_size)
            : m_nx(nx)
            , m_ny(ny)
            , m_nz(nz)
            , m_cell_size(cell_size)
            , m_velocity(nx * ny * nz, Vec3{})
            , m_pressure(nx * ny * nz, 0.0)
            , m_divergence(nx * ny * nz, 0.0)
            , m_liquid_fraction(nx * ny * nz, 1.0)
            , m_gas_fraction(nx * ny * nz, 0.0)
            , m_curvature(nx * ny * nz, 0.0)
            , m_gradient(nx * ny * nz, Vec3{})
            , m_weights(nx * ny * nz, 0.0)
        {
        }

        [[nodiscard]] std::size_t nx() const noexcept { return m_nx; }
        [[nodiscard]] std::size_t ny() const noexcept { return m_ny; }
        [[nodiscard]] std::size_t nz() const noexcept { return m_nz; }
        [[nodiscard]] double cell_size() const noexcept { return m_cell_size; }

        [[nodiscard]] std::size_t index(std::size_t ix, std::size_t iy, std::size_t iz) const noexcept
        {
            return ix + m_nx * (iy + m_ny * iz);
        }

        [[nodiscard]] Vec3 cell_center(std::size_t ix, std::size_t iy, std::size_t iz) const noexcept
        {
            return Vec3{
                (static_cast<double>(ix) + 0.5) * m_cell_size,
                (static_cast<double>(iy) + 0.5) * m_cell_size,
                (static_cast<double>(iz) + 0.5) * m_cell_size,
            };
        }

        void clear()
        {
            std::fill(m_velocity.begin(), m_velocity.end(), Vec3{});
            std::fill(m_pressure.begin(), m_pressure.end(), 0.0);
            std::fill(m_divergence.begin(), m_divergence.end(), 0.0);
            std::fill(m_liquid_fraction.begin(), m_liquid_fraction.end(), 1.0);
            std::fill(m_gas_fraction.begin(), m_gas_fraction.end(), 0.0);
            std::fill(m_curvature.begin(), m_curvature.end(), 0.0);
            std::fill(m_gradient.begin(), m_gradient.end(), Vec3{});
            std::fill(m_weights.begin(), m_weights.end(), 0.0);
        }

        [[nodiscard]] Vec3 sample_velocity(const Vec3& position) const noexcept
        {
            return sample_vector_field(m_velocity, position);
        }

        [[nodiscard]] Vec3 sample_gradient(const Vec3& position) const noexcept
        {
            return sample_vector_field(m_gradient, position);
        }

        [[nodiscard]] double sample_pressure(const Vec3& position) const noexcept
        {
            return sample_scalar_field(m_pressure, position);
        }

        [[nodiscard]] double sample_gas_fraction(const Vec3& position) const noexcept
        {
            return sample_scalar_field(m_gas_fraction, position);
        }

        [[nodiscard]] double sample_curvature(const Vec3& position) const noexcept
        {
            return sample_scalar_field(m_curvature, position);
        }

        std::vector<Vec3>& velocities() noexcept { return m_velocity; }
        const std::vector<Vec3>& velocities() const noexcept { return m_velocity; }

        std::vector<double>& pressure() noexcept { return m_pressure; }
        const std::vector<double>& pressure() const noexcept { return m_pressure; }

        std::vector<double>& divergence() noexcept { return m_divergence; }
        const std::vector<double>& divergence() const noexcept { return m_divergence; }

        std::vector<double>& liquid_fraction() noexcept { return m_liquid_fraction; }
        const std::vector<double>& liquid_fraction() const noexcept { return m_liquid_fraction; }

        std::vector<double>& gas_fraction() noexcept { return m_gas_fraction; }
        const std::vector<double>& gas_fraction() const noexcept { return m_gas_fraction; }

        std::vector<double>& curvature() noexcept { return m_curvature; }
        const std::vector<double>& curvature() const noexcept { return m_curvature; }

        std::vector<Vec3>& gradient() noexcept { return m_gradient; }
        const std::vector<Vec3>& gradient() const noexcept { return m_gradient; }

        std::vector<double>& weights() noexcept { return m_weights; }
        const std::vector<double>& weights() const noexcept { return m_weights; }

    private:
        [[nodiscard]] Vec3 sample_vector_field(const std::vector<Vec3>& field, const Vec3& position) const noexcept
        {
            const double fx = std::clamp(position.x / m_cell_size - 0.5, 0.0, static_cast<double>(m_nx - 1));
            const double fy = std::clamp(position.y / m_cell_size - 0.5, 0.0, static_cast<double>(m_ny - 1));
            const double fz = std::clamp(position.z / m_cell_size - 0.5, 0.0, static_cast<double>(m_nz - 1));

            const std::size_t ix0 = static_cast<std::size_t>(std::floor(fx));
            const std::size_t iy0 = static_cast<std::size_t>(std::floor(fy));
            const std::size_t iz0 = static_cast<std::size_t>(std::floor(fz));

            const std::size_t ix1 = std::min(ix0 + 1, m_nx - 1);
            const std::size_t iy1 = std::min(iy0 + 1, m_ny - 1);
            const std::size_t iz1 = std::min(iz0 + 1, m_nz - 1);

            const double tx = std::clamp(fx - static_cast<double>(ix0), 0.0, 1.0);
            const double ty = std::clamp(fy - static_cast<double>(iy0), 0.0, 1.0);
            const double tz = std::clamp(fz - static_cast<double>(iz0), 0.0, 1.0);

            const auto idx000 = index(ix0, iy0, iz0);
            const auto idx100 = index(ix1, iy0, iz0);
            const auto idx010 = index(ix0, iy1, iz0);
            const auto idx110 = index(ix1, iy1, iz0);
            const auto idx001 = index(ix0, iy0, iz1);
            const auto idx101 = index(ix1, iy0, iz1);
            const auto idx011 = index(ix0, iy1, iz1);
            const auto idx111 = index(ix1, iy1, iz1);

            const Vec3 c00 = field[idx000] * (1.0 - tx) + field[idx100] * tx;
            const Vec3 c10 = field[idx010] * (1.0 - tx) + field[idx110] * tx;
            const Vec3 c01 = field[idx001] * (1.0 - tx) + field[idx101] * tx;
            const Vec3 c11 = field[idx011] * (1.0 - tx) + field[idx111] * tx;

            const Vec3 c0 = c00 * (1.0 - ty) + c10 * ty;
            const Vec3 c1 = c01 * (1.0 - ty) + c11 * ty;

            return c0 * (1.0 - tz) + c1 * tz;
        }

        [[nodiscard]] double sample_scalar_field(const std::vector<double>& field, const Vec3& position) const noexcept
        {
            const double fx = std::clamp(position.x / m_cell_size - 0.5, 0.0, static_cast<double>(m_nx - 1));
            const double fy = std::clamp(position.y / m_cell_size - 0.5, 0.0, static_cast<double>(m_ny - 1));
            const double fz = std::clamp(position.z / m_cell_size - 0.5, 0.0, static_cast<double>(m_nz - 1));

            const std::size_t ix0 = static_cast<std::size_t>(std::floor(fx));
            const std::size_t iy0 = static_cast<std::size_t>(std::floor(fy));
            const std::size_t iz0 = static_cast<std::size_t>(std::floor(fz));

            const std::size_t ix1 = std::min(ix0 + 1, m_nx - 1);
            const std::size_t iy1 = std::min(iy0 + 1, m_ny - 1);
            const std::size_t iz1 = std::min(iz0 + 1, m_nz - 1);

            const double tx = std::clamp(fx - static_cast<double>(ix0), 0.0, 1.0);
            const double ty = std::clamp(fy - static_cast<double>(iy0), 0.0, 1.0);
            const double tz = std::clamp(fz - static_cast<double>(iz0), 0.0, 1.0);

            const auto idx000 = index(ix0, iy0, iz0);
            const auto idx100 = index(ix1, iy0, iz0);
            const auto idx010 = index(ix0, iy1, iz0);
            const auto idx110 = index(ix1, iy1, iz0);
            const auto idx001 = index(ix0, iy0, iz1);
            const auto idx101 = index(ix1, iy0, iz1);
            const auto idx011 = index(ix0, iy1, iz1);
            const auto idx111 = index(ix1, iy1, iz1);

            const double c00 = field[idx000] * (1.0 - tx) + field[idx100] * tx;
            const double c10 = field[idx010] * (1.0 - tx) + field[idx110] * tx;
            const double c01 = field[idx001] * (1.0 - tx) + field[idx101] * tx;
            const double c11 = field[idx011] * (1.0 - tx) + field[idx111] * tx;

            const double c0 = c00 * (1.0 - ty) + c10 * ty;
            const double c1 = c01 * (1.0 - ty) + c11 * ty;

            return c0 * (1.0 - tz) + c1 * tz;
        }

        std::size_t m_nx{0};
        std::size_t m_ny{0};
        std::size_t m_nz{0};
        double m_cell_size{1.0};

        std::vector<Vec3> m_velocity{};
        std::vector<double> m_pressure{};
        std::vector<double> m_divergence{};
        std::vector<double> m_liquid_fraction{};
        std::vector<double> m_gas_fraction{};
        std::vector<double> m_curvature{};
        std::vector<Vec3> m_gradient{};
        std::vector<double> m_weights{};
    };

    class BubbleSimulation
    {
    public:
        BubbleSimulation(std::size_t nx, std::size_t ny, std::size_t nz, double cell_size)
            : m_grid(nx, ny, nz, cell_size)
        {
        }

        void add_bubble(const Vec3& position, double radius)
        {
            Bubble bubble{};
            bubble.position = position;
            bubble.radius = radius;
            bubble.mass = std::max(1e-6, m_bubble_density * bubble.volume());
            m_bubbles.push_back(bubble);
        }

        [[nodiscard]] const std::vector<Bubble>& bubbles() const noexcept { return m_bubbles; }
        [[nodiscard]] std::vector<Bubble>& bubbles() noexcept { return m_bubbles; }

        [[nodiscard]] const EulerianGrid& grid() const noexcept { return m_grid; }
        [[nodiscard]] EulerianGrid& grid() noexcept { return m_grid; }

        void step(double dt)
        {
            if (m_bubbles.empty() || dt <= std::numeric_limits<double>::epsilon())
            {
                return;
            }

            m_grid.clear();

            deposit_bubble_velocities();
            average_grid_velocities();
            compute_volume_fractions();
            compute_surface_tension_fields();
            compute_divergence();
            solve_pressure(dt);
            apply_pressure_gradient(dt);
            apply_forces(dt);
            resolve_collisions();
            clamp_bubbles_to_domain();
        }

    private:
        void deposit_bubble_velocities()
        {
            auto& velocities = m_grid.velocities();
            auto& weights = m_grid.weights();

            for (const auto& bubble : m_bubbles)
            {
                const double inv_cell = 1.0 / m_grid.cell_size();
                const int ix = static_cast<int>(std::clamp(bubble.position.x * inv_cell, 0.0, static_cast<double>(m_grid.nx() - 1)));
                const int iy = static_cast<int>(std::clamp(bubble.position.y * inv_cell, 0.0, static_cast<double>(m_grid.ny() - 1)));
                const int iz = static_cast<int>(std::clamp(bubble.position.z * inv_cell, 0.0, static_cast<double>(m_grid.nz() - 1)));

                const auto idx = m_grid.index(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy), static_cast<std::size_t>(iz));
                velocities[idx] += bubble.velocity * bubble.volume();
                weights[idx] += bubble.volume();
            }
        }

        void average_grid_velocities()
        {
            auto& velocities = m_grid.velocities();
            auto& weights = m_grid.weights();

            for (std::size_t i = 0; i < velocities.size(); ++i)
            {
                if (weights[i] > std::numeric_limits<double>::epsilon())
                {
                    velocities[i] /= weights[i];
                }
            }
        }
        }

        void compute_volume_fractions()
        {
            auto& gas_fraction = m_grid.gas_fraction();
            auto& liquid_fraction = m_grid.liquid_fraction();
            std::fill(gas_fraction.begin(), gas_fraction.end(), 0.0);
            std::fill(liquid_fraction.begin(), liquid_fraction.end(), 1.0);

            const double cell_volume = std::pow(m_grid.cell_size(), 3.0);

            for (const auto& bubble : m_bubbles)
            {
                const double radius = bubble.radius;
                const int min_ix = std::max(0, static_cast<int>(std::floor((bubble.position.x - radius) / m_grid.cell_size())));
                const int min_iy = std::max(0, static_cast<int>(std::floor((bubble.position.y - radius) / m_grid.cell_size())));
                const int min_iz = std::max(0, static_cast<int>(std::floor((bubble.position.z - radius) / m_grid.cell_size())));
                const int max_ix = std::min(static_cast<int>(m_grid.nx()) - 1, static_cast<int>(std::ceil((bubble.position.x + radius) / m_grid.cell_size())));
                const int max_iy = std::min(static_cast<int>(m_grid.ny()) - 1, static_cast<int>(std::ceil((bubble.position.y + radius) / m_grid.cell_size())));
                const int max_iz = std::min(static_cast<int>(m_grid.nz()) - 1, static_cast<int>(std::ceil((bubble.position.z + radius) / m_grid.cell_size())));

                double total_weight = 0.0;
                std::vector<std::pair<std::size_t, double>> contributions;
                contributions.reserve(static_cast<std::size_t>((max_ix - min_ix + 1) * (max_iy - min_iy + 1) * (max_iz - min_iz + 1)));

                for (int iz = min_iz; iz <= max_iz; ++iz)
                {
                    for (int iy = min_iy; iy <= max_iy; ++iy)
                    {
                        for (int ix = min_ix; ix <= max_ix; ++ix)
                        {
                            const auto center = m_grid.cell_center(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy), static_cast<std::size_t>(iz));
                            const Vec3 delta{center.x - bubble.position.x, center.y - bubble.position.y, center.z - bubble.position.z};
                            const double distance = length(delta);
                            if (distance > radius)
                            {
                                continue;
                            }
                            const double weight = 1.0 - (distance / std::max(radius, 1e-6));
                            if (weight <= 0.0)
                            {
                                continue;
                            }
                            const auto idx = m_grid.index(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy), static_cast<std::size_t>(iz));
                            contributions.emplace_back(idx, weight);
                            total_weight += weight;
                        }
                    }
                }

                if (contributions.empty() || total_weight <= std::numeric_limits<double>::epsilon())
                {
                    continue;
                }

                const double scale = std::min(1.0, bubble.volume() / (cell_volume * total_weight));
                for (const auto& [index, weight] : contributions)
                {
                    gas_fraction[index] += scale * weight;
                }
            }

            for (std::size_t i = 0; i < gas_fraction.size(); ++i)
            {
                gas_fraction[i] = std::clamp(gas_fraction[i], 0.0, 1.0);
                liquid_fraction[i] = std::clamp(1.0 - gas_fraction[i], 0.0, 1.0);
            }
        }

        void compute_surface_tension_fields()
        {
            auto& gradient = m_grid.gradient();
            auto& curvature = m_grid.curvature();
            const auto& gas_fraction = m_grid.gas_fraction();

            const double inv_cell = 1.0 / m_grid.cell_size();
            const std::size_t nx = m_grid.nx();
            const std::size_t ny = m_grid.ny();
            const std::size_t nz = m_grid.nz();

            for (std::size_t k = 0; k < nz; ++k)
            {
                for (std::size_t j = 0; j < ny; ++j)
                {
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);

                        const auto sample = [&](int x, int y, int z) {
                            x = std::clamp(x, 0, static_cast<int>(nx) - 1);
                            y = std::clamp(y, 0, static_cast<int>(ny) - 1);
                            z = std::clamp(z, 0, static_cast<int>(nz) - 1);
                            return gas_fraction[m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z))];
                        };

                        const double gx = (sample(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k)) - sample(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k))) * 0.5 * inv_cell;
                        const double gy = (sample(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k)) - sample(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k))) * 0.5 * inv_cell;
                        const double gz = (sample(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1) - sample(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1)) * 0.5 * inv_cell;

                        gradient[idx] = Vec3{gx, gy, gz};
                    }
                }
            }

            for (std::size_t k = 0; k < nz; ++k)
            {
                for (std::size_t j = 0; j < ny; ++j)
                {
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);

                        const auto sample_normal = [&](int x, int y, int z) {
                            x = std::clamp(x, 0, static_cast<int>(nx) - 1);
                            y = std::clamp(y, 0, static_cast<int>(ny) - 1);
                            z = std::clamp(z, 0, static_cast<int>(nz) - 1);
                            const auto id = m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z));
                            return normalise(gradient[id]);
                        };

                        const Vec3 nxp = sample_normal(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 nxm = sample_normal(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 nyp = sample_normal(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        const Vec3 nym = sample_normal(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        const Vec3 nzp = sample_normal(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        const Vec3 nzm = sample_normal(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);

                        const double div_n = (nxp.x - nxm.x + nyp.y - nym.y + nzp.z - nzm.z) * 0.5 * inv_cell;
                        curvature[idx] = -div_n;
                    }
                }
            }
        }

        void compute_divergence()
        {
            auto& divergence = m_grid.divergence();
            const auto& velocity = m_grid.velocities();

            const double inv_cell = 1.0 / m_grid.cell_size();
            const std::size_t nx = m_grid.nx();
            const std::size_t ny = m_grid.ny();
            const std::size_t nz = m_grid.nz();

            for (std::size_t k = 0; k < nz; ++k)
            {
                for (std::size_t j = 0; j < ny; ++j)
                {
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);

                        const auto sample_velocity = [&](int x, int y, int z) {
                            x = std::clamp(x, 0, static_cast<int>(nx) - 1);
                            y = std::clamp(y, 0, static_cast<int>(ny) - 1);
                            z = std::clamp(z, 0, static_cast<int>(nz) - 1);
                            return velocity[m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z))];
                        };

                        const Vec3 vx_p = sample_velocity(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 vx_m = sample_velocity(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 vy_p = sample_velocity(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        const Vec3 vy_m = sample_velocity(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        const Vec3 vz_p = sample_velocity(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        const Vec3 vz_m = sample_velocity(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);

                        const double dudx = (vx_p.x - vx_m.x) * 0.5 * inv_cell;
                        const double dvdy = (vy_p.y - vy_m.y) * 0.5 * inv_cell;
                        const double dwdz = (vz_p.z - vz_m.z) * 0.5 * inv_cell;

                        divergence[idx] = (dudx + dvdy + dwdz) * m_grid.liquid_fraction()[idx];
                    }
                }
            }
        }
        void solve_pressure(double dt)
        {
            const std::size_t nx = m_grid.nx();
            const std::size_t ny = m_grid.ny();
            const std::size_t nz = m_grid.nz();
            const std::size_t cell_count = nx * ny * nz;

            std::vector<std::vector<std::size_t>> adjacency(cell_count);
            std::vector<bool> is_dirichlet(cell_count, false);

            const auto& liquid_fraction = m_grid.liquid_fraction();
            const double threshold = 0.05;

            for (std::size_t k = 0; k < nz; ++k)
            {
                for (std::size_t j = 0; j < ny; ++j)
                {
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);
                        adjacency[idx].push_back(idx);

                        if (liquid_fraction[idx] <= threshold)
                        {
                            is_dirichlet[idx] = true;
                            continue;
                        }

                        const auto add_neighbor = [&](int x, int y, int z) {
                            if (x < 0 || y < 0 || z < 0 || x >= static_cast<int>(nx) || y >= static_cast<int>(ny) || z >= static_cast<int>(nz))
                            {
                                return;
                            }
                            const auto neighbor_idx = m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z));
                            adjacency[idx].push_back(neighbor_idx);
                        };

                        add_neighbor(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        add_neighbor(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        add_neighbor(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        add_neighbor(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        add_neighbor(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        add_neighbor(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);
                    }
                }
            }

            for (std::size_t row = 0; row < cell_count; ++row)
            {
                auto& entries = adjacency[row];
                std::sort(entries.begin(), entries.end());
                entries.erase(std::unique(entries.begin(), entries.end()), entries.end());
            }

            std::vector<int> row_prefix(cell_count + 1, 0);
            for (std::size_t row = 0; row < cell_count; ++row)
            {
                row_prefix[row + 1] = row_prefix[row] + static_cast<int>(adjacency[row].size());
            }

            detail::CooMatrix coo{};
            coo.dimension = cell_count;
            coo.row_prefix = row_prefix;
            coo.rows.assign(static_cast<std::size_t>(coo.row_prefix.back()), 0);
            coo.cols.assign(static_cast<std::size_t>(coo.row_prefix.back()), 0);
            coo.values.assign(static_cast<std::size_t>(coo.row_prefix.back()), 0.0);

            std::vector<std::unordered_map<std::size_t, std::size_t>> index_map(cell_count);
            for (std::size_t row = 0; row < cell_count; ++row)
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
                    index_map[row].emplace(column, index);
                }
            }

            std::vector<double> rhs(cell_count, 0.0);
            const auto& divergence = m_grid.divergence();
            const double inv_h2 = 1.0 / (m_grid.cell_size() * m_grid.cell_size());

            for (std::size_t k = 0; k < nz; ++k)
            {
                for (std::size_t j = 0; j < ny; ++j)
                {
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto row = m_grid.index(i, j, k);

                        if (is_dirichlet[row])
                        {
                            const auto diag_it = index_map[row].find(row);
                            if (diag_it == index_map[row].end())
                            {
                                throw std::logic_error("Missing diagonal for Dirichlet cell");
                            }
                            coo.values[diag_it->second] = 1.0;
                            rhs[row] = 0.0;
                            continue;
                        }

                        double diagonal = 0.0;

                        const auto handle_neighbor = [&](int x, int y, int z) {
                            if (x < 0 || y < 0 || z < 0 || x >= static_cast<int>(nx) || y >= static_cast<int>(ny) || z >= static_cast<int>(nz))
                            {
                                diagonal += inv_h2;
                                return;
                            }

                            const auto neighbor = m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z));
                            const double weight = 0.5 * (liquid_fraction[row] + liquid_fraction[neighbor]);
                            if (weight <= threshold)
                            {
                                diagonal += inv_h2;
                                return;
                            }

                            const auto entry = index_map[row].find(neighbor);
                            if (entry != index_map[row].end())
                            {
                                coo.values[entry->second] -= inv_h2;
                            }

                            diagonal += inv_h2;
                        };

                        handle_neighbor(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        handle_neighbor(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        handle_neighbor(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        handle_neighbor(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        handle_neighbor(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        handle_neighbor(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);

                        const auto diag_it = index_map[row].find(row);
                        if (diag_it == index_map[row].end())
                        {
                            throw std::logic_error("Missing diagonal entry");
                        }

                        coo.values[diag_it->second] += diagonal;
                        rhs[row] = (m_fluid_density / dt) * divergence[row];
                    }
                }
            }

            detail::CsrMatrix csr(coo);

            SolverOptions options{};
            options.solver = SolverType::ConjugateGradient;
            options.preconditioner = PreconditionerType::Jacobi;
            options.tolerance = 1e-6;
            options.max_iterations = static_cast<std::size_t>(std::max<std::size_t>(1000, cell_count * 4));

            const auto summary = detail::solve_linear_system(csr, rhs, options);
            m_grid.pressure() = summary.solution;
        }
        void apply_pressure_gradient(double dt)
        {
            auto& velocity = m_grid.velocities();
            const auto& pressure = m_grid.pressure();

            const double inv_density = 1.0 / m_fluid_density;
            const double inv_cell = 1.0 / m_grid.cell_size();
            const std::size_t nx = m_grid.nx();
            const std::size_t ny = m_grid.ny();
            const std::size_t nz = m_grid.nz();

            for (std::size_t k = 0; k < nz; ++k)
            {
                for (std::size_t j = 0; j < ny; ++j)
                {
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);

                        const auto sample_pressure = [&](int x, int y, int z) {
                            x = std::clamp(x, 0, static_cast<int>(nx) - 1);
                            y = std::clamp(y, 0, static_cast<int>(ny) - 1);
                            z = std::clamp(z, 0, static_cast<int>(nz) - 1);
                            return pressure[m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z))];
                        };

                        const double dpdx = (sample_pressure(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k)) - sample_pressure(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k))) * 0.5 * inv_cell;
                        const double dpdy = (sample_pressure(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k)) - sample_pressure(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k))) * 0.5 * inv_cell;
                        const double dpdz = (sample_pressure(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1) - sample_pressure(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1)) * 0.5 * inv_cell;

                        velocity[idx].x -= dt * inv_density * dpdx;
                        velocity[idx].y -= dt * inv_density * dpdy;
                        velocity[idx].z -= dt * inv_density * dpdz;
                    }
                }
            }
        }

        void apply_forces(double dt)
        {
            for (auto& bubble : m_bubbles)
            {
                const Vec3 fluid_velocity = m_grid.sample_velocity(bubble.position);
                const Vec3 relative_velocity = bubble.velocity - fluid_velocity;
                const Vec3 drag_force = -m_drag_coefficient * relative_velocity * bubble.surface_area();

                const double gas_fraction = m_grid.sample_gas_fraction(bubble.position);
                const double buoyancy_strength = (1.0 - gas_fraction) * bubble.volume() * (m_fluid_density - m_bubble_density) * m_gravity;
                const Vec3 buoyancy_force{0.0, buoyancy_strength, 0.0};

                const Vec3 gradient = m_grid.sample_gradient(bubble.position);
                Vec3 normal = normalise(gradient);
                const double curvature = m_grid.sample_curvature(bubble.position);
                Vec3 surface_force{};
                if (length(normal) > 0.0)
                {
                    surface_force = -m_surface_tension * curvature * normal * bubble.surface_area();
                }

                const Vec3 total_force = drag_force + buoyancy_force + surface_force;
                const Vec3 acceleration = total_force / bubble.mass;
                bubble.velocity += acceleration * dt;

                bubble.velocity = bubble.velocity * (1.0 - m_coupling_strength) + fluid_velocity * m_coupling_strength;
                bubble.position += bubble.velocity * dt;
            }
        }

        void resolve_collisions()
        {
            if (m_bubbles.size() < 2)
            {
                return;
            }

            const double restitution = 0.5;

            for (std::size_t i = 0; i < m_bubbles.size(); ++i)
            {
                for (std::size_t j = i + 1; j < m_bubbles.size(); ++j)
                {
                    auto& a = m_bubbles[i];
                    auto& b = m_bubbles[j];

                    const Vec3 delta{b.position.x - a.position.x, b.position.y - a.position.y, b.position.z - a.position.z};
                    const double dist = length(delta);
                    const double min_dist = a.radius + b.radius;

                    if (dist < min_dist && dist > std::numeric_limits<double>::epsilon())
                    {
                        const Vec3 normal = delta / dist;
                        const double penetration = min_dist - dist;

                        a.position -= normal * (0.5 * penetration);
                        b.position += normal * (0.5 * penetration);

                        const Vec3 relative_velocity{b.velocity.x - a.velocity.x, b.velocity.y - a.velocity.y, b.velocity.z - a.velocity.z};
                        const double rel_normal = dot(relative_velocity, normal);
                        if (rel_normal > 0.0)
                        {
                            continue;
                        }

                        const double impulse = -(1.0 + restitution) * rel_normal / (1.0 / a.mass + 1.0 / b.mass);
                        const Vec3 impulse_vec = impulse * normal;
                        a.velocity -= impulse_vec / a.mass;
                        b.velocity += impulse_vec / b.mass;
                    }
                }
            }
        }

        void clamp_bubbles_to_domain()
        {
            const double max_x = static_cast<double>(m_grid.nx()) * m_grid.cell_size();
            const double max_y = static_cast<double>(m_grid.ny()) * m_grid.cell_size();
            const double max_z = static_cast<double>(m_grid.nz()) * m_grid.cell_size();

            for (auto& bubble : m_bubbles)
            {
                const auto clamp_axis = [&](double& value, double& velocity, double radius, double max_bound) {
                    if (value < radius)
                    {
                        value = radius;
                        velocity *= -0.5;
                    }
                    else if (value > max_bound - radius)
                    {
                        value = max_bound - radius;
                        velocity *= -0.5;
                    }
                };

                clamp_axis(bubble.position.x, bubble.velocity.x, bubble.radius, max_x);
                clamp_axis(bubble.position.y, bubble.velocity.y, bubble.radius, max_y);
                clamp_axis(bubble.position.z, bubble.velocity.z, bubble.radius, max_z);
            }
        }

        EulerianGrid m_grid;
        std::vector<Bubble> m_bubbles{};

        double m_fluid_density{1000.0};
        double m_bubble_density{1.2};
        double m_gravity{9.81};
        double m_drag_coefficient{6.0};
        double m_surface_tension{0.0728};
        double m_coupling_strength{0.2};
    };
} // namespace almond::fem::bubbles

