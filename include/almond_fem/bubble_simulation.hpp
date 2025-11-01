#pragma once

// Harden against Windows.h min/max macros seen elsewhere in the TU.
#if defined(_WIN32) && !defined(NOMINMAX)
#define NOMINMAX
#endif

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
        double x{ 0.0 };
        double y{ 0.0 };
        double z{ 0.0 };

        Vec3& operator+=(const Vec3& other) noexcept { x += other.x; y += other.y; z += other.z; return *this; }
        Vec3& operator-=(const Vec3& other) noexcept { x -= other.x; y -= other.y; z -= other.z; return *this; }
        Vec3& operator*=(double s) noexcept { x *= s; y *= s; z *= s; return *this; }
        Vec3& operator/=(double s) {
            if (std::abs(s) <= std::numeric_limits<double>::epsilon()) throw std::runtime_error("Vec3 / 0");
            x /= s; y /= s; z /= s; return *this;
        }
    };

    [[nodiscard]] inline Vec3 operator+(Vec3 a, const Vec3& b) noexcept { a += b; return a; }
    [[nodiscard]] inline Vec3 operator-(Vec3 a, const Vec3& b) noexcept { a -= b; return a; }
    [[nodiscard]] inline Vec3 operator*(Vec3 a, double s) noexcept { a *= s; return a; }
    [[nodiscard]] inline Vec3 operator*(double s, Vec3 b) noexcept { b *= s; return b; }
    [[nodiscard]] inline Vec3 operator/(Vec3 a, double s) { a /= s; return a; }

    [[nodiscard]] inline double dot(const Vec3& a, const Vec3& b) noexcept { return a.x * b.x + a.y * b.y + a.z * b.z; }
    [[nodiscard]] inline double length_squared(const Vec3& v) noexcept { return dot(v, v); }
    [[nodiscard]] inline double length(const Vec3& v) noexcept { return std::sqrt(length_squared(v)); }

    [[nodiscard]] inline Vec3 normalise(const Vec3& v)
    {
        const double len = length(v);
        if (len <= std::numeric_limits<double>::epsilon()) return {};
        return v / len;
    }

    struct Bubble
    {
        Vec3 position{};
        Vec3 velocity{};
        double radius{ 0.5 };
        double mass{ 1.0 };

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
            : m_nx(nx), m_ny(ny), m_nz(nz), m_cell_size(cell_size),
            m_velocity(nx* ny* nz, Vec3{}),
            m_pressure(nx* ny* nz, 0.0),
            m_divergence(nx* ny* nz, 0.0),
            m_liquid_fraction(nx* ny* nz, 1.0),
            m_gas_fraction(nx* ny* nz, 0.0),
            m_curvature(nx* ny* nz, 0.0),
            m_gradient(nx* ny* nz, Vec3{}),
            m_weights(nx* ny* nz, 0.0)
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
            return {
                (static_cast<double>(ix) + 0.5) * m_cell_size,
                (static_cast<double>(iy) + 0.5) * m_cell_size,
                (static_cast<double>(iz) + 0.5) * m_cell_size
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

        [[nodiscard]] Vec3 sample_velocity(const Vec3& p) const noexcept { return sample_vector_field(m_velocity, p); }
        [[nodiscard]] Vec3 sample_gradient(const Vec3& p) const noexcept { return sample_vector_field(m_gradient, p); }
        [[nodiscard]] double sample_pressure(const Vec3& p) const noexcept { return sample_scalar_field(m_pressure, p); }
        [[nodiscard]] double sample_gas_fraction(const Vec3& p) const noexcept { return sample_scalar_field(m_gas_fraction, p); }
        [[nodiscard]] double sample_curvature(const Vec3& p) const noexcept { return sample_scalar_field(m_curvature, p); }

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
        [[nodiscard]] Vec3 sample_vector_field(const std::vector<Vec3>& field, const Vec3& pos) const noexcept
        {
            const double fx = (std::clamp)(pos.x / m_cell_size - 0.5, 0.0, static_cast<double>(m_nx - 1));
            const double fy = (std::clamp)(pos.y / m_cell_size - 0.5, 0.0, static_cast<double>(m_ny - 1));
            const double fz = (std::clamp)(pos.z / m_cell_size - 0.5, 0.0, static_cast<double>(m_nz - 1));

            const std::size_t ix0 = static_cast<std::size_t>(std::floor(fx));
            const std::size_t iy0 = static_cast<std::size_t>(std::floor(fy));
            const std::size_t iz0 = static_cast<std::size_t>(std::floor(fz));

            const std::size_t ix1 = (std::min)(ix0 + 1, m_nx - 1);
            const std::size_t iy1 = (std::min)(iy0 + 1, m_ny - 1);
            const std::size_t iz1 = (std::min)(iz0 + 1, m_nz - 1);

            const double tx = (std::clamp)(fx - static_cast<double>(ix0), 0.0, 1.0);
            const double ty = (std::clamp)(fy - static_cast<double>(iy0), 0.0, 1.0);
            const double tz = (std::clamp)(fz - static_cast<double>(iz0), 0.0, 1.0);

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

        [[nodiscard]] double sample_scalar_field(const std::vector<double>& field, const Vec3& pos) const noexcept
        {
            const double fx = (std::clamp)(pos.x / m_cell_size - 0.5, 0.0, static_cast<double>(m_nx - 1));
            const double fy = (std::clamp)(pos.y / m_cell_size - 0.5, 0.0, static_cast<double>(m_ny - 1));
            const double fz = (std::clamp)(pos.z / m_cell_size - 0.5, 0.0, static_cast<double>(m_nz - 1));

            const std::size_t ix0 = static_cast<std::size_t>(std::floor(fx));
            const std::size_t iy0 = static_cast<std::size_t>(std::floor(fy));
            const std::size_t iz0 = static_cast<std::size_t>(std::floor(fz));

            const std::size_t ix1 = (std::min)(ix0 + 1, m_nx - 1);
            const std::size_t iy1 = (std::min)(iy0 + 1, m_ny - 1);
            const std::size_t iz1 = (std::min)(iz0 + 1, m_nz - 1);

            const double tx = (std::clamp)(fx - static_cast<double>(ix0), 0.0, 1.0);
            const double ty = (std::clamp)(fy - static_cast<double>(iy0), 0.0, 1.0);
            const double tz = (std::clamp)(fz - static_cast<double>(iz0), 0.0, 1.0);

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

        std::size_t m_nx{ 0 }, m_ny{ 0 }, m_nz{ 0 };
        double m_cell_size{ 1.0 };

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
            : m_grid(nx, ny, nz, cell_size) {
        }

        void add_bubble(const Vec3& position, double radius)
        {
            Bubble b{};
            b.position = position;
            b.radius = radius;
            // Avoid Windows macro traps with (std::max)
            b.mass = (std::max)(1e-6, m_bubble_density * b.volume());
            m_bubbles.push_back(b);
        }

        [[nodiscard]] const std::vector<Bubble>& bubbles() const noexcept { return m_bubbles; }
        [[nodiscard]] std::vector<Bubble>& bubbles() noexcept { return m_bubbles; }

        [[nodiscard]] const EulerianGrid& grid() const noexcept { return m_grid; }
        [[nodiscard]] EulerianGrid& grid() noexcept { return m_grid; }

        void step(double dt)
        {
            if (m_bubbles.empty() || dt <= std::numeric_limits<double>::epsilon()) return;

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
            auto& vel = m_grid.velocities();
            auto& w = m_grid.weights();

            for (const auto& b : m_bubbles)
            {
                const double inv_cell = 1.0 / m_grid.cell_size();
                const int ix = static_cast<int>((std::clamp)(b.position.x * inv_cell, 0.0, static_cast<double>(m_grid.nx() - 1)));
                const int iy = static_cast<int>((std::clamp)(b.position.y * inv_cell, 0.0, static_cast<double>(m_grid.ny() - 1)));
                const int iz = static_cast<int>((std::clamp)(b.position.z * inv_cell, 0.0, static_cast<double>(m_grid.nz() - 1)));

                const auto idx = m_grid.index(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy), static_cast<std::size_t>(iz));
                vel[idx] += b.velocity * b.volume();
                w[idx] += b.volume();
            }
        }

        void average_grid_velocities()
        {
            auto& vel = m_grid.velocities();
            auto& w = m_grid.weights();
            for (std::size_t i = 0; i < vel.size(); ++i)
                if (w[i] > std::numeric_limits<double>::epsilon())
                    vel[i] /= w[i];
        }

        void compute_volume_fractions()
        {
            auto& gas = m_grid.gas_fraction();
            auto& liq = m_grid.liquid_fraction();
            std::fill(gas.begin(), gas.end(), 0.0);
            std::fill(liq.begin(), liq.end(), 1.0);

            const double h = m_grid.cell_size();
            const double cell_volume = h * h * h;

            for (const auto& b : m_bubbles)
            {
                const double r = b.radius;

                const int min_ix = (std::max)(0, static_cast<int>(std::floor((b.position.x - r) / h)));
                const int min_iy = (std::max)(0, static_cast<int>(std::floor((b.position.y - r) / h)));
                const int min_iz = (std::max)(0, static_cast<int>(std::floor((b.position.z - r) / h)));

                const int max_ix = (std::min)(static_cast<int>(m_grid.nx()) - 1, static_cast<int>(std::ceil((b.position.x + r) / h)));
                const int max_iy = (std::min)(static_cast<int>(m_grid.ny()) - 1, static_cast<int>(std::ceil((b.position.y + r) / h)));
                const int max_iz = (std::min)(static_cast<int>(m_grid.nz()) - 1, static_cast<int>(std::ceil((b.position.z + r) / h)));

                double total_w = 0.0;
                std::vector<std::pair<std::size_t, double>> contrib;
                contrib.reserve(static_cast<std::size_t>((max_ix - min_ix + 1) * (max_iy - min_iy + 1) * (max_iz - min_iz + 1)));

                for (int iz = min_iz; iz <= max_iz; ++iz)
                    for (int iy = min_iy; iy <= max_iy; ++iy)
                        for (int ix = min_ix; ix <= max_ix; ++ix)
                        {
                            const auto c = m_grid.cell_center(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy), static_cast<std::size_t>(iz));
                            const Vec3 d{ c.x - b.position.x, c.y - b.position.y, c.z - b.position.z };
                            const double dist = length(d);
                            if (dist > r) continue;

                            const double w = 1.0 - (dist / (std::max)(r, 1e-6));
                            if (w <= 0.0) continue;

                            const auto idx = m_grid.index(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy), static_cast<std::size_t>(iz));
                            contrib.emplace_back(idx, w);
                            total_w += w;
                        }

                if (contrib.empty() || total_w <= std::numeric_limits<double>::epsilon()) continue;

                const double scale = (std::min)(1.0, b.volume() / (cell_volume * total_w));
                for (const auto& [idx, w] : contrib) gas[idx] += scale * w;
            }

            for (std::size_t i = 0; i < gas.size(); ++i)
            {
                gas[i] = (std::clamp)(gas[i], 0.0, 1.0);
                liq[i] = (std::clamp)(1.0 - gas[i], 0.0, 1.0);
            }
        }

        void compute_surface_tension_fields()
        {
            auto& grad = m_grid.gradient();
            auto& curv = m_grid.curvature();
            const auto& gas = m_grid.gas_fraction();

            const double inv_h = 1.0 / m_grid.cell_size();
            const std::size_t nx = m_grid.nx(), ny = m_grid.ny(), nz = m_grid.nz();

            // Gradient of gas fraction
            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t j = 0; j < ny; ++j)
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);
                        const auto sample = [&](int x, int y, int z) {
                            x = (std::clamp)(x, 0, static_cast<int>(nx) - 1);
                            y = (std::clamp)(y, 0, static_cast<int>(ny) - 1);
                            z = (std::clamp)(z, 0, static_cast<int>(nz) - 1);
                            return gas[m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z))];
                            };
                        const double gx = (sample(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k)) -
                            sample(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k))) * 0.5 * inv_h;
                        const double gy = (sample(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k)) -
                            sample(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k))) * 0.5 * inv_h;
                        const double gz = (sample(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1) -
                            sample(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1)) * 0.5 * inv_h;
                        grad[idx] = Vec3{ gx, gy, gz };
                    }

            // Curvature = -div(n), n = normalize(grad)
            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t j = 0; j < ny; ++j)
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);
                        const auto n_at = [&](int x, int y, int z) {
                            x = (std::clamp)(x, 0, static_cast<int>(nx) - 1);
                            y = (std::clamp)(y, 0, static_cast<int>(ny) - 1);
                            z = (std::clamp)(z, 0, static_cast<int>(nz) - 1);
                            return normalise(m_grid.gradient()[m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z))]);
                            };
                        const Vec3 nxp = n_at(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 nxm = n_at(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 nyp = n_at(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        const Vec3 nym = n_at(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        const Vec3 nzp = n_at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        const Vec3 nzm = n_at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);

                        const double div_n = (nxp.x - nxm.x + nyp.y - nym.y + nzp.z - nzm.z) * 0.5 * inv_h;
                        curv[idx] = -div_n;
                    }
        }

        void compute_divergence()
        {
            auto& div = m_grid.divergence();
            const auto& vel = m_grid.velocities();

            const double inv_h = 1.0 / m_grid.cell_size();
            const std::size_t nx = m_grid.nx(), ny = m_grid.ny(), nz = m_grid.nz();

            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t j = 0; j < ny; ++j)
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);
                        const auto v_at = [&](int x, int y, int z) {
                            x = (std::clamp)(x, 0, static_cast<int>(nx) - 1);
                            y = (std::clamp)(y, 0, static_cast<int>(ny) - 1);
                            z = (std::clamp)(z, 0, static_cast<int>(nz) - 1);
                            return vel[m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z))];
                            };
                        const Vec3 vxp = v_at(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 vxm = v_at(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        const Vec3 vyp = v_at(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        const Vec3 vym = v_at(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        const Vec3 vzp = v_at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        const Vec3 vzm = v_at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);

                        const double dudx = (vxp.x - vxm.x) * 0.5 * inv_h;
                        const double dvdy = (vyp.y - vym.y) * 0.5 * inv_h;
                        const double dwdz = (vzp.z - vzm.z) * 0.5 * inv_h;

                        div[idx] = (dudx + dvdy + dwdz) * m_grid.liquid_fraction()[idx];
                    }
        }

        void solve_pressure(double dt)
        {
            const std::size_t nx = m_grid.nx(), ny = m_grid.ny(), nz = m_grid.nz();
            const std::size_t N = nx * ny * nz;

            std::vector<std::vector<std::size_t>> adj(N);
            std::vector<bool> is_dirichlet(N, false);

            const auto& liq = m_grid.liquid_fraction();
            const double threshold = 0.05;

            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t j = 0; j < ny; ++j)
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);
                        adj[idx].push_back(idx);

                        if (liq[idx] <= threshold) { is_dirichlet[idx] = true; continue; }

                        const auto push_n = [&](int x, int y, int z) {
                            if (x < 0 || y < 0 || z < 0 ||
                                x >= static_cast<int>(nx) ||
                                y >= static_cast<int>(ny) ||
                                z >= static_cast<int>(nz)) return;
                            adj[idx].push_back(m_grid.index(static_cast<std::size_t>(x),
                                static_cast<std::size_t>(y),
                                static_cast<std::size_t>(z)));
                            };
                        push_n(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        push_n(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        push_n(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        push_n(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        push_n(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        push_n(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);
                    }

            for (auto& row : adj) { std::sort(row.begin(), row.end()); row.erase(std::unique(row.begin(), row.end()), row.end()); }

            std::vector<int> row_prefix(N + 1, 0);
            for (std::size_t r = 0; r < N; ++r) row_prefix[r + 1] = row_prefix[r] + static_cast<int>(adj[r].size());

            detail::CooMatrix coo{};
            coo.dimension = N;
            coo.row_prefix = row_prefix;
            coo.rows.assign(static_cast<std::size_t>(coo.row_prefix.back()), 0);
            coo.cols.assign(static_cast<std::size_t>(coo.row_prefix.back()), 0);
            coo.values.assign(static_cast<std::size_t>(coo.row_prefix.back()), 0.0);

            std::vector<std::unordered_map<std::size_t, std::size_t>> map(N);
            for (std::size_t r = 0; r < N; ++r)
            {
                const auto start = static_cast<std::size_t>(coo.row_prefix[r]);
                const auto& cols = adj[r];
                map[r].reserve(cols.size());
                for (std::size_t off = 0; off < cols.size(); ++off)
                {
                    const auto idx = start + off;
                    const auto c = cols[off];
                    coo.rows[idx] = static_cast<int>(r);
                    coo.cols[idx] = static_cast<int>(c);
                    map[r].emplace(c, idx);
                }
            }

            std::vector<double> rhs(N, 0.0);
            const auto& div = m_grid.divergence();
            const double inv_h2 = 1.0 / (m_grid.cell_size() * m_grid.cell_size());

            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t j = 0; j < ny; ++j)
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto row = m_grid.index(i, j, k);

                        if (is_dirichlet[row])
                        {
                            const auto it = map[row].find(row);
                            if (it == map[row].end()) throw std::logic_error("Missing diagonal for Dirichlet cell");
                            coo.values[it->second] = 1.0;
                            rhs[row] = 0.0;
                            continue;
                        }

                        double diag = 0.0;
                        const auto handle = [&](int x, int y, int z) {
                            if (x < 0 || y < 0 || z < 0 ||
                                x >= static_cast<int>(nx) ||
                                y >= static_cast<int>(ny) ||
                                z >= static_cast<int>(nz))
                            {
                                diag += inv_h2; return;
                            }

                            const auto n = m_grid.index(static_cast<std::size_t>(x),
                                static_cast<std::size_t>(y),
                                static_cast<std::size_t>(z));
                            const double w = 0.5 * (liq[row] + liq[n]);
                            if (w <= threshold) { diag += inv_h2; return; }

                            const auto e = map[row].find(n);
                            if (e != map[row].end()) coo.values[e->second] -= inv_h2;
                            diag += inv_h2;
                            };

                        handle(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k));
                        handle(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k));
                        handle(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k));
                        handle(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k));
                        handle(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1);
                        handle(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1);

                        const auto it = map[row].find(row);
                        if (it == map[row].end()) throw std::logic_error("Missing diagonal entry");

                        coo.values[it->second] += diag;
                        rhs[row] = (m_fluid_density / dt) * div[row];
                    }

            detail::CsrMatrix csr(coo);

            SolverOptions opt{};
            opt.solver = SolverType::ConjugateGradient;
            opt.preconditioner = PreconditionerType::Jacobi;
            opt.tolerance = 1e-6;
            opt.max_iterations = static_cast<std::size_t>((std::max<std::size_t>)(1000, N * 4));

            const auto summary = detail::solve_linear_system(csr, rhs, opt);
            m_grid.pressure() = summary.solution;
        }

        void apply_pressure_gradient(double dt)
        {
            auto& vel = m_grid.velocities();
            const auto& p = m_grid.pressure();

            const double inv_rho = 1.0 / m_fluid_density;
            const double inv_h = 1.0 / m_grid.cell_size();
            const std::size_t nx = m_grid.nx(), ny = m_grid.ny(), nz = m_grid.nz();

            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t j = 0; j < ny; ++j)
                    for (std::size_t i = 0; i < nx; ++i)
                    {
                        const auto idx = m_grid.index(i, j, k);
                        const auto p_at = [&](int x, int y, int z) {
                            x = (std::clamp)(x, 0, static_cast<int>(nx) - 1);
                            y = (std::clamp)(y, 0, static_cast<int>(ny) - 1);
                            z = (std::clamp)(z, 0, static_cast<int>(nz) - 1);
                            return p[m_grid.index(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<std::size_t>(z))];
                            };
                        const double dpdx = (p_at(static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k)) -
                            p_at(static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k))) * 0.5 * inv_h;
                        const double dpdy = (p_at(static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k)) -
                            p_at(static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k))) * 0.5 * inv_h;
                        const double dpdz = (p_at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1) -
                            p_at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1)) * 0.5 * inv_h;

                        vel[idx].x -= dt * inv_rho * dpdx;
                        vel[idx].y -= dt * inv_rho * dpdy;
                        vel[idx].z -= dt * inv_rho * dpdz;
                    }
        }

        void apply_forces(double dt)
        {
            for (auto& b : m_bubbles)
            {
                const Vec3 u = m_grid.sample_velocity(b.position);
                const Vec3 rel = b.velocity - u;
                const Vec3 Fdrag = -m_drag_coefficient * rel * b.surface_area();

                const double gas = m_grid.sample_gas_fraction(b.position);
                const double buoy = (1.0 - gas) * b.volume() * (m_fluid_density - m_bubble_density) * m_gravity;
                const Vec3 Fbuoy{ 0.0, buoy, 0.0 };

                const Vec3 grad = m_grid.sample_gradient(b.position);
                const Vec3 n = normalise(grad);
                const double kappa = m_grid.sample_curvature(b.position);
                Vec3 Fsurf{};
                if (length(n) > 0.0) Fsurf = -m_surface_tension * kappa * n * b.surface_area();

                const Vec3 F = Fdrag + Fbuoy + Fsurf;
                const Vec3 a = F / b.mass;
                b.velocity += a * dt;

                // Coupling
                b.velocity = b.velocity * (1.0 - m_coupling_strength) + u * m_coupling_strength;
                b.position += b.velocity * dt;
            }
        }

        void resolve_collisions()
        {
            if (m_bubbles.size() < 2) return;

            const double restitution = 0.5;
            for (std::size_t i = 0; i < m_bubbles.size(); ++i)
                for (std::size_t j = i + 1; j < m_bubbles.size(); ++j)
                {
                    auto& a = m_bubbles[i];
                    auto& b = m_bubbles[j];

                    const Vec3 d{ b.position.x - a.position.x, b.position.y - a.position.y, b.position.z - a.position.z };
                    const double dist = length(d);
                    const double min_dist = a.radius + b.radius;

                    if (dist < min_dist && dist > std::numeric_limits<double>::epsilon())
                    {
                        const Vec3 n = d / dist;
                        const double pen = min_dist - dist;

                        a.position -= n * (0.5 * pen);
                        b.position += n * (0.5 * pen);

                        const Vec3 rel{ b.velocity.x - a.velocity.x, b.velocity.y - a.velocity.y, b.velocity.z - a.velocity.z };
                        const double rel_n = dot(rel, n);
                        if (rel_n > 0.0) continue;

                        const double J = -(1.0 + restitution) * rel_n / (1.0 / a.mass + 1.0 / b.mass);
                        const Vec3 Jn = J * n;
                        a.velocity -= Jn / a.mass;
                        b.velocity += Jn / b.mass;
                    }
                }
        }

        void clamp_bubbles_to_domain()
        {
            const double max_x = static_cast<double>(m_grid.nx()) * m_grid.cell_size();
            const double max_y = static_cast<double>(m_grid.ny()) * m_grid.cell_size();
            const double max_z = static_cast<double>(m_grid.nz()) * m_grid.cell_size();

            for (auto& b : m_bubbles)
            {
                const auto clamp_axis = [](double& v, double& vv, double r, double maxb) {
                    if (v < r) { v = r; vv *= -0.5; }
                    else if (v > maxb - r) { v = maxb - r; vv *= -0.5; }
                    };
                clamp_axis(b.position.x, b.velocity.x, b.radius, max_x);
                clamp_axis(b.position.y, b.velocity.y, b.radius, max_y);
                clamp_axis(b.position.z, b.velocity.z, b.radius, max_z);
            }
        }

    private:
        EulerianGrid m_grid;
        std::vector<Bubble> m_bubbles{};

        double m_fluid_density{ 1000.0 };
        double m_bubble_density{ 1.2 };
        double m_gravity{ 9.81 };
        double m_drag_coefficient{ 6.0 };
        double m_surface_tension{ 0.0728 };
        double m_coupling_strength{ 0.2 };
    };
} // namespace almond::fem::bubbles
