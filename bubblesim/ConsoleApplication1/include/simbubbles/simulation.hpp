#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>

// Basic 3D vector struct for convenience
struct Vec3 {
    float x, y, z;
    Vec3(float x_ = 0, float y_ = 0, float z_ = 0) : x(x_), y(y_), z(z_) {}

    Vec3 operator-() const { return Vec3(-x, -y, -z); }

    Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(float s) const { return Vec3(x / s, y / s, z / s); }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator*=(float s) { x *= s; y *= s; z *= s; return *this; }

    float dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    float norm() const { return std::sqrt(x * x + y * y + z * z); }
    Vec3 normalized() const { float n = norm(); return n > 1e-8f ? (*this) * (1.0f / n) : Vec3(); }
};

// Particle representing a bubble
struct Particle {
    Vec3 pos;
    Vec3 vel;
    float radius;
    float volume;
    float mass;
    Particle(const Vec3& p, float r, float rho_b)
        : pos(p), vel(0, 0, 0), radius(r) {
        volume = (4.0f / 3.0f) * 3.1415926f * r * r * r;
        mass = rho_b * volume;
    }
};

class FluidSim {
public:
    // Grid dimensions
    int Nx, Ny, Nz;
    float dx;               // cell size (assuming cubic cells)
    // Fluid properties
    float rho_w;            // water density
    float rho_b;            // bubble (gas) density
    float mu;               // dynamic viscosity of water
    float sigma;            // surface tension coefficient
    Vec3 g;                 // gravity acceleration (m/s^2, negative in y or z direction)
    // Simulation state
    std::vector<Particle> particles;
    // MAC grid storage (staggered):
    // u-component velocity at faces between i and i+1 in x-direction
    std::vector<float> u;
    // v-component velocity at faces between j and j+1 in y-direction
    std::vector<float> v;
    // w-component velocity at faces between k and k+1 in z-direction
    std::vector<float> w;
    // Pressure at cell centers
    std::vector<float> pressure;
    // Bubble volume fraction at cells
    std::vector<float> phi_b;
    // Temporary array for previous phi_b (to compute changes)
    std::vector<float> phi_b_prev;
    // Divergence (for pressure solve)
    std::vector<float> divergence;

    FluidSim(int nx, int ny, int nz, float cell_size)
        : Nx(nx), Ny(ny), Nz(nz), dx(cell_size) {
        // Initialize fluid properties
        rho_w = 1000.0f;    // water density (kg/m^3)
        rho_b = 1.0f;       // bubble gas density (kg/m^3)
        mu = 0.001f;     // water viscosity (Pa·s, ~water)
        sigma = 0.072f;     // surface tension (N/m, ~water-air)
        // Gravity vector (pointing downward in Y direction here)
        g = Vec3(0.0f, -9.81f, 0.0f);
        // Allocate grid arrays with proper sizing
        u.resize((Nx + 1) * Ny * Nz);
        v.resize(Nx * (Ny + 1) * Nz);
        w.resize(Nx * Ny * (Nz + 1));
        pressure.resize(Nx * Ny * Nz);
        phi_b.resize(Nx * Ny * Nz);
        phi_b_prev.resize(Nx * Ny * Nz);
        divergence.resize(Nx * Ny * Nz);
        // Initialize fields to zero
        std::fill(u.begin(), u.end(), 0.0f);
        std::fill(v.begin(), v.end(), 0.0f);
        std::fill(w.begin(), w.end(), 0.0f);
        std::fill(pressure.begin(), pressure.end(), 0.0f);
        std::fill(phi_b.begin(), phi_b.end(), 0.0f);
        std::fill(phi_b_prev.begin(), phi_b_prev.end(), 0.0f);
    }

    // Utility: index calculations for flattened arrays
    inline int idxCell(int i, int j, int k) const {
        return (k * Ny + j) * Nx + i;
    }
    inline int idxU(int i, int j, int k) const {
        // u at face (i between cells i-1 and i, j, k)
        return (k * Ny + j) * (Nx + 1) + i;
    }
    inline int idxV(int i, int j, int k) const {
        // v at face (i, j between j-1 and j, k)
        return (k * (Ny + 1) + j) * Nx + i;
    }
    inline int idxW(int i, int j, int k) const {
        // w at face (i, j, k between k-1 and k)
        return (k * Ny + j) * Nx + i;
    }

    // Add a bubble particle to the simulation
    void addBubble(const Vec3& position, float radius) {
        particles.emplace_back(position, radius, rho_b);
    }

    // Interpolate water velocity at an arbitrary point (within domain)
    Vec3 sampleVelocity(const Vec3& pos) const {
        float fx = pos.x / dx;
        float fy = pos.y / dx;
        float fz = pos.z / dx;

        int i0 = static_cast<int>(floor(fx));
        int j0 = static_cast<int>(floor(fy));
        int k0 = static_cast<int>(floor(fz));

        float fracx = fx - i0;
        float fracy = fy - j0;
        float fracz = fz - k0;

        // Clamp indices
        i0 = std::clamp(i0, 0, Nx - 1);
        j0 = std::clamp(j0, 0, Ny - 1);
        k0 = std::clamp(k0, 0, Nz - 1);

        // --- U-component (x-faces) ---
        int iU = std::clamp(i0, 0, Nx);
        int iU1 = std::clamp(iU + 1, 0, Nx);
        int jU = j0;
        int kU = k0;

        float ux000 = u[idxU(iU, jU, kU)];
        float ux010 = u[idxU(iU, jU + 1 < Ny ? jU + 1 : Ny - 1, kU)];
        float ux001 = u[idxU(iU, jU, kU + 1 < Nz ? kU + 1 : Nz - 1)];
        float ux011 = u[idxU(iU, jU + 1 < Ny ? jU + 1 : Ny - 1, kU + 1 < Nz ? kU + 1 : Nz - 1)];
        float ux100 = u[idxU(iU1, jU, kU)];
        float ux110 = u[idxU(iU1, jU + 1 < Ny ? jU + 1 : Ny - 1, kU)];
        float ux101 = u[idxU(iU1, jU, kU + 1 < Nz ? kU + 1 : Nz - 1)];
        float ux111 = u[idxU(iU1, jU + 1 < Ny ? jU + 1 : Ny - 1, kU + 1 < Nz ? kU + 1 : Nz - 1)];

        float ux00 = lerp(ux000, ux010, fracy);
        float ux01 = lerp(ux001, ux011, fracy);
        float ux10 = lerp(ux100, ux110, fracy);
        float ux11 = lerp(ux101, ux111, fracy);
        float ux0 = lerp(ux00, ux01, fracz);
        float ux1 = lerp(ux10, ux11, fracz);
        float u_val = lerp(ux0, ux1, fracx);

        // --- V-component (y-faces) ---
        int jV = std::clamp(j0, 0, Ny);
        int jV1 = std::clamp(jV + 1, 0, Ny);
        int iV = i0;
        int kV = k0;

        float vy000 = v[idxV(iV, jV, kV)];
        float vy010 = v[idxV(iV, jV1, kV)];
        float vy001 = v[idxV(iV, jV, kV + 1 < Nz ? kV + 1 : Nz - 1)];
        float vy011 = v[idxV(iV, jV1, kV + 1 < Nz ? kV + 1 : Nz - 1)];
        float vy100 = v[idxV(iV + 1 < Nx ? iV + 1 : Nx - 1, jV, kV)];
        float vy110 = v[idxV(iV + 1 < Nx ? iV + 1 : Nx - 1, jV1, kV)];
        float vy101 = v[idxV(iV + 1 < Nx ? iV + 1 : Nx - 1, jV, kV + 1 < Nz ? kV + 1 : Nz - 1)];
        float vy111 = v[idxV(iV + 1 < Nx ? iV + 1 : Nx - 1, jV1, kV + 1 < Nz ? kV + 1 : Nz - 1)];

        float vy00 = lerp(vy000, vy010, fracy);
        float vy01 = lerp(vy001, vy011, fracy);
        float vy10 = lerp(vy100, vy110, fracy);
        float vy11 = lerp(vy101, vy111, fracy);
        float vy0 = lerp(vy00, vy01, fracz);
        float vy1 = lerp(vy10, vy11, fracz);
        float v_val = lerp(vy0, vy1, fracx);

        // --- W-component (z-faces) ---
        int kW = std::clamp(k0, 0, Nz);
        int kW1 = std::clamp(kW + 1, 0, Nz);
        int iW = i0;
        int jW = j0;

        float wz000 = w[idxW(iW, jW, kW)];
        float wz010 = w[idxW(iW, jW + 1 < Ny ? jW + 1 : Ny - 1, kW)];
        float wz001 = w[idxW(iW, jW, kW1)];
        float wz011 = w[idxW(iW, jW + 1 < Ny ? jW + 1 : Ny - 1, kW1)];
        float wz100 = w[idxW(iW + 1 < Nx ? iW + 1 : Nx - 1, jW, kW)];
        float wz110 = w[idxW(iW + 1 < Nx ? iW + 1 : Nx - 1, jW + 1 < Ny ? jW + 1 : Ny - 1, kW)];
        float wz101 = w[idxW(iW + 1 < Nx ? iW + 1 : Nx - 1, jW, kW1)];
        float wz111 = w[idxW(iW + 1 < Nx ? iW + 1 : Nx - 1, jW + 1 < Ny ? jW + 1 : Ny - 1, kW1)];

        float wz00 = lerp(wz000, wz010, fracy);
        float wz01 = lerp(wz001, wz011, fracy);
        float wz10 = lerp(wz100, wz110, fracy);
        float wz11 = lerp(wz101, wz111, fracy);
        float wz0 = lerp(wz00, wz01, fracz);
        float wz1 = lerp(wz10, wz11, fracz);
        float w_val = lerp(wz0, wz1, fracx);

        return Vec3(u_val, v_val, w_val);
    }

    // Deposit a force (per unit volume) into the fluid velocity field (affecting water momentum)
    void depositForce(const Vec3& pos, const Vec3& f) {
        float fx = pos.x / dx;
        float fy = pos.y / dx;
        float fz = pos.z / dx;

        int i0 = static_cast<int>(floor(fx));
        int j0 = static_cast<int>(floor(fy));
        int k0 = static_cast<int>(floor(fz));

        float fracx = fx - i0;
        float fracy = fy - j0;
        float fracz = fz - k0;

        // Clamp indices
        i0 = std::clamp(i0, 0, Nx - 1);
        j0 = std::clamp(j0, 0, Ny - 1);
        k0 = std::clamp(k0, 0, Nz - 1);

        float rho = rho_w;

        // Precompute weights
        float w000 = (1 - fracx) * (1 - fracy) * (1 - fracz);
        float w100 = fracx * (1 - fracy) * (1 - fracz);
        float w010 = (1 - fracx) * fracy * (1 - fracz);
        float w110 = fracx * fracy * (1 - fracz);
        float w001 = (1 - fracx) * (1 - fracy) * fracz;
        float w101 = fracx * (1 - fracy) * fracz;
        float w011 = (1 - fracx) * fracy * fracz;
        float w111 = fracx * fracy * fracz;

        // ----- U component (x faces) -----
        int iU = i0;
        int jU = j0;
        int kU = k0;

        u[idxU(iU, jU, kU)] += (f.x / rho) * w000;
        if (jU + 1 < Ny) u[idxU(iU, jU + 1, kU)] += (f.x / rho) * w010;
        if (kU + 1 < Nz) u[idxU(iU, jU, kU + 1)] += (f.x / rho) * w001;
        if (jU + 1 < Ny && kU + 1 < Nz) u[idxU(iU, jU + 1, kU + 1)] += (f.x / rho) * w011;

        if (iU + 1 < Nx + 1) {
            u[idxU(iU + 1, jU, kU)] += (f.x / rho) * w100;
            if (jU + 1 < Ny) u[idxU(iU + 1, jU + 1, kU)] += (f.x / rho) * w110;
            if (kU + 1 < Nz) u[idxU(iU + 1, jU, kU + 1)] += (f.x / rho) * w101;
            if (jU + 1 < Ny && kU + 1 < Nz) u[idxU(iU + 1, jU + 1, kU + 1)] += (f.x / rho) * w111;
        }

        // ----- V component (y faces) -----
        int iV = i0;
        int jV = j0;
        int kV = k0;

        v[idxV(iV, jV, kV)] += (f.y / rho) * w000;
        if (iV + 1 < Nx) v[idxV(iV + 1, jV, kV)] += (f.y / rho) * w100;
        if (kV + 1 < Nz) v[idxV(iV, jV, kV + 1)] += (f.y / rho) * w001;
        if (iV + 1 < Nx && kV + 1 < Nz) v[idxV(iV + 1, jV, kV + 1)] += (f.y / rho) * w101;

        if (jV + 1 < Ny + 1) {
            v[idxV(iV, jV + 1, kV)] += (f.y / rho) * w010;
            if (iV + 1 < Nx) v[idxV(iV + 1, jV + 1, kV)] += (f.y / rho) * w110;
            if (kV + 1 < Nz) v[idxV(iV, jV + 1, kV + 1)] += (f.y / rho) * w011;
            if (iV + 1 < Nx && kV + 1 < Nz) v[idxV(iV + 1, jV + 1, kV + 1)] += (f.y / rho) * w111;
        }

        // ----- W component (z faces) -----
        int iW = i0;
        int jW = j0;
        int kW = k0;

        w[idxW(iW, jW, kW)] += (f.z / rho) * w000;
        if (iW + 1 < Nx) w[idxW(iW + 1, jW, kW)] += (f.z / rho) * w100;
        if (jW + 1 < Ny) w[idxW(iW, jW + 1, kW)] += (f.z / rho) * w010;
        if (iW + 1 < Nx && jW + 1 < Ny) w[idxW(iW + 1, jW + 1, kW)] += (f.z / rho) * w110;

        if (kW + 1 < Nz + 1) {
            w[idxW(iW, jW, kW + 1)] += (f.z / rho) * w001;
            if (iW + 1 < Nx) w[idxW(iW + 1, jW, kW + 1)] += (f.z / rho) * w101;
            if (jW + 1 < Ny) w[idxW(iW, jW + 1, kW + 1)] += (f.z / rho) * w011;
            if (iW + 1 < Nx && jW + 1 < Ny) w[idxW(iW + 1, jW + 1, kW + 1)] += (f.z / rho) * w111;
        }
    }


    // Advect bubble particles (simple explicit Euler integration)
    void advectBubbles(float dt) {
        for (auto& p : particles) {
            // Update position: x_new = x_old + v * dt
            p.pos += p.vel * dt;
            // Keep particle within bounds (reflect or clamp at walls)
            // Here we implement simple elastic reflection at domain boundaries
            if (p.pos.x < 0) { p.pos.x = 0; if (p.vel.x < 0) p.vel.x *= -0.5f; } // damp a bit to simulate energy loss
            if (p.pos.x > Nx * dx) { p.pos.x = Nx * dx; if (p.vel.x > 0) p.vel.x *= -0.5f; }
            if (p.pos.y < 0) { p.pos.y = 0; if (p.vel.y < 0) p.vel.y *= -0.5f; }
            if (p.pos.y > Ny * dx) { p.pos.y = Ny * dx; if (p.vel.y > 0) p.vel.y *= -0.5f; }
            if (p.pos.z < 0) { p.pos.z = 0; if (p.vel.z < 0) p.vel.z *= -0.5f; }
            if (p.pos.z > Nz * dx) { p.pos.z = Nz * dx; if (p.vel.z > 0) p.vel.z *= -0.5f; }
        }
    }

    // Update φ_b field from particle volumes
    void updatePhiB() {
        // Save previous φ_b
        phi_b_prev = phi_b;
        // Reset current φ_b
        std::fill(phi_b.begin(), phi_b.end(), 0.0f);
        // Deposit each particle's volume into its containing cell
        for (auto& p : particles) {
            int i = static_cast<int>(p.pos.x / dx);
            int j = static_cast<int>(p.pos.y / dx);
            int k = static_cast<int>(p.pos.z / dx);
            if (i < 0) i = 0;
            if (i >= Nx) i = Nx - 1;
            if (j < 0) j = 0;
            if (j >= Ny) j = Ny - 1;
            if (k < 0) k = 0;
            if (k >= Nz) k = Nz - 1;
            int ci = idxCell(i, j, k);
            // Add bubble volume fraction = particle volume / cell volume
            float cellVol = dx * dx * dx;
            float volFraction = p.volume / cellVol;
            phi_b[ci] += volFraction;
            if (phi_b[ci] > 1.0f) {
                phi_b[ci] = 1.0f; // clamp if overfull
            }
        }
        // Compute water fraction implicitly as 1-phi_b when needed (not stored explicitly)
    }

    // Compute surface tension force per face using φ_b field (Brackbill's CSF approach)
    void computeSurfaceTension(std::vector<Vec3>& f_sigma_faces) {
        f_sigma_faces.clear();
        f_sigma_faces.resize((Nx + 1) * Ny * Nz + Nx * (Ny + 1) * Nz + Nx * Ny * (Nz + 1)); // allocate for all faces (we won't use all in one vector, but for simplicity)
        // Compute normals and curvature at cell centers
        std::vector<Vec3> normals;
        normals.resize(Nx * Ny * Nz);
        std::vector<float> curvature;
        curvature.resize(Nx * Ny * Nz);
        // Compute normalized gradient of φ_b at cell centers
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    // central differences (Neumann boundary: use one-sided at edges)
                    float phi_c = phi_b[idxCell(i, j, k)];
                    float phi_x0 = (i > 0 ? phi_b[idxCell(i - 1, j, k)] : phi_c);
                    float phi_x1 = (i < Nx - 1 ? phi_b[idxCell(i + 1, j, k)] : phi_c);
                    float phi_y0 = (j > 0 ? phi_b[idxCell(i, j - 1, k)] : phi_c);
                    float phi_y1 = (j < Ny - 1 ? phi_b[idxCell(i, j + 1, k)] : phi_c);
                    float phi_z0 = (k > 0 ? phi_b[idxCell(i, j, k - 1)] : phi_c);
                    float phi_z1 = (k < Nz - 1 ? phi_b[idxCell(i, j, k + 1)] : phi_c);
                    Vec3 grad;
                    grad.x = (phi_x1 - phi_x0) / dx;
                    grad.y = (phi_y1 - phi_y0) / dx;
                    grad.z = (phi_z1 - phi_z0) / dx;
                    normals[idxCell(i, j, k)] = grad;
                }
            }
        }
        // Smooth normals (optional small smoothing to reduce noise)
        for (int idx = 0; idx < (int)normals.size(); ++idx) {
            normals[idx] = normals[idx].normalized();
        }
        // Compute curvature = divergence of normals
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    Vec3 n_c = normals[idxCell(i, j, k)];
                    // divergence: take neighbor normals difference
                    float nx0 = (i > 0 ? normals[idxCell(i - 1, j, k)].x : n_c.x);
                    float nx1 = (i < Nx - 1 ? normals[idxCell(i + 1, j, k)].x : n_c.x);
                    float ny0 = (j > 0 ? normals[idxCell(i, j - 1, k)].y : n_c.y);
                    float ny1 = (j < Ny - 1 ? normals[idxCell(i, j + 1, k)].y : n_c.y);
                    float nz0 = (k > 0 ? normals[idxCell(i, j, k - 1)].z : n_c.z);
                    float nz1 = (k < Nz - 1 ? normals[idxCell(i, j, k + 1)].z : n_c.z);
                    float div_n = (nx1 - nx0) / (2 * dx) + (ny1 - ny0) / (2 * dx) + (nz1 - nz0) / (2 * dx);
                    curvature[idxCell(i, j, k)] = div_n;
                }
            }
        }
        // Compute face-centered surface tension forces
        // X-faces
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i <= Nx; ++i) {
                    // faces in x-direction exist between cell (i-1) and (i)
                    if (i == 0 || i == Nx) {
                        // at domain boundary, no surface tension force (assuming no air outside or treat as symmetry)
                        f_sigma_faces[idxU(i, j, k)] = Vec3(0, 0, 0);
                    }
                    else {
                        int cellL = idxCell(i - 1, j, k);
                        int cellR = idxCell(i, j, k);
                        // difference in φ_b across face
                        float dphi = phi_b[cellR] - phi_b[cellL];
                        // average curvature at face
                        float kappa_avg = 0.5f * (curvature[cellL] + curvature[cellR]);
                        float f_scalar = sigma * kappa_avg * dphi / dx;
                        // Force direction: along +x (from left to right) when φ_b difference is positive
                        Vec3 f_face(f_scalar, 0.0f, 0.0f);
                        // assign force per volume for this face (acting on fluid)
                        f_sigma_faces[idxU(i, j, k)] = f_face;
                    }
                }
            }
        }
        // Y-faces
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j <= Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    if (j == 0 || j == Ny) {
                        f_sigma_faces[idxV(i, j, k)] = Vec3(0, 0, 0);
                    }
                    else {
                        int cellB = idxCell(i, j - 1, k);
                        int cellT = idxCell(i, j, k);
                        float dphi = phi_b[cellT] - phi_b[cellB];
                        float kappa_avg = 0.5f * (curvature[cellB] + curvature[cellT]);
                        float f_scalar = sigma * kappa_avg * dphi / dx;
                        // direction +y
                        Vec3 f_face(0.0f, f_scalar, 0.0f);
                        f_sigma_faces[idxV(i, j, k)] = f_face;
                    }
                }
            }
        }
        // Z-faces
        for (int k = 0; k <= Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    if (k == 0 || k == Nz) {
                        f_sigma_faces[idxW(i, j, k)] = Vec3(0, 0, 0);
                    }
                    else {
                        int cellF = idxCell(i, j, k - 1);
                        int cellB = idxCell(i, j, k);
                        float dphi = phi_b[cellB] - phi_b[cellF];
                        float kappa_avg = 0.5f * (curvature[cellF] + curvature[cellB]);
                        float f_scalar = sigma * kappa_avg * dphi / dx;
                        // direction +z
                        Vec3 f_face(0.0f, 0.0f, f_scalar);
                        f_sigma_faces[idxW(i, j, k)] = f_face;
                    }
                }
            }
        }
    }

    // Apply forces: gravity, buoyancy, drag, surface tension
    void applyForces(float dt) {
        // Gravity & Buoyancy:
        // Water: apply gravity to water velocity (will be balanced by pressure in incompressible solve)
        // Here we choose not to directly accelerate water with gravity to avoid accumulating drift; 
        // Instead, gravity's effect on water will appear via buoyancy pressure differences (hydrostatic).
        // Bubbles: apply gravity (downward) and buoyancy (upward)
        for (auto& p : particles) {
            // Gravity acceleration on bubble
            Vec3 a_grav = g; // g is downward (e.g. (0,-9.81,0))
            // Buoyancy acceleration = (rho_w - rho_b)/rho_b * g (upward relative to gravity)
            float addedMassFactor = 0.5f; // added mass of water (approx 0.5 of displaced mass for a sphere)
            float effectiveRho = rho_b + addedMassFactor * rho_w;
            Vec3 a_buoy = g * (-rho_w / effectiveRho); // negative of gravity * (rho_w/eff_rho) gives upward accel
            // Net acceleration from gravity+buoyancy:
            Vec3 a_net = a_grav + a_buoy;
            p.vel += a_net * dt;
        }

        // Drag:
        for (auto& p : particles) {
            // Interpolate water velocity at particle position
            Vec3 u_w = sampleVelocity(p.pos);
            Vec3 rel = u_w - p.vel; // relative velocity (water minus bubble)
            // Only compute drag if there's significant relative motion or water present
            if (rel.norm() > 1e-6f && phi_w_at(p.pos) > 1e-3f) {
                // Compute linear + quadratic drag
                float r = p.radius;
                // Linear Stokes drag: F_lin = 6 π μ r * rel
                Vec3 F_lin = rel * (6.0f * 3.1415926f * mu * r);
                // Quadratic drag: F_quad = 0.5 * C_d * ρ_w * A * ||rel|| * rel  (A = π r^2)
                float C_d = 0.5f; // drag coefficient
                float A = 3.1415926f * r * r;
                Vec3 F_quad = rel.normalized() * (0.5f * C_d * rho_w * A * rel.norm() * rel.norm());
                Vec3 F_drag = F_lin + F_quad;
                // Scale drag by water fraction at particle (no drag in pure air or pure water edge cases)
                float phi_w_part = phi_w_at(p.pos);
                F_drag *= phi_w_part;
                // Apply drag to bubble (particle)
                Vec3 a_drag = F_drag / p.mass;
                p.vel += a_drag * dt;
                // Apply opposite drag to water (momentum transfer)
                Vec3 f_drag = -F_drag / (dx * dx * dx); // force per unit volume (approx, using cell volume)
                depositForce(p.pos, f_drag * 1.0f); // deposit into fluid velocities
            }
        }

        // Surface tension:
        // Compute face forces due to surface tension
        std::vector<Vec3> f_sigma;
        computeSurfaceTension(f_sigma);
        // Apply a fraction of surface tension force to avoid instability (explicit integration)
        float STF = 0.2f; // surface tension force scaling (less than 1 for stability)
        // Add to water velocity field
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 1; i < Nx; ++i) { // x faces (skip boundaries)
                    Vec3 f = f_sigma[idxU(i, j, k)] * STF;
                    // water impulse: Δu = (f/ρ_w)*dt (only water fraction gets this force)
                    u[idxU(i, j, k)] += (f.x / rho_w) * dt;
                }
            }
        }
        for (int k = 0; k < Nz; ++k) {
            for (int j = 1; j < Ny; ++j) { // y faces
                for (int i = 0; i < Nx; ++i) {
                    Vec3 f = f_sigma[idxV(i, j, k)] * STF;
                    v[idxV(i, j, k)] += (f.y / rho_w) * dt;
                }
            }
        }
        for (int k = 1; k < Nz; ++k) { // z faces
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    Vec3 f = f_sigma[idxW(i, j, k)] * STF;
                    w[idxW(i, j, k)] += (f.z / rho_w) * dt;
                }
            }
        }
    }

    // Helper: compute water fraction at a particle's position (by trilinear interpolation of φ_w)
    float phi_w_at(const Vec3& pos) const {
        // Compute cell indices and fractions similar to velocity sampling
        float fx = pos.x / dx;
        float fy = pos.y / dx;
        float fz = pos.z / dx;
        int i0 = static_cast<int>(floor(fx));
        int j0 = static_cast<int>(floor(fy));
        int k0 = static_cast<int>(floor(fz));
        float fracx = fx - i0;
        float fracy = fy - j0;
        float fracz = fz - k0;
        if (i0 < 0) { i0 = 0; fracx = 0; }
        if (j0 < 0) { j0 = 0; fracy = 0; }
        if (k0 < 0) { k0 = 0; fracz = 0; }
        if (i0 > Nx - 2) { i0 = Nx - 2; fracx = 1; }
        if (j0 > Ny - 2) { j0 = Ny - 2; fracy = 1; }
        if (k0 > Nz - 2) { k0 = Nz - 2; fracz = 1; }
        // gather φ_b at eight surrounding cells
        float c000 = phi_b[idxCell(i0, j0, k0)];
        float c100 = phi_b[idxCell(i0 + 1, j0, k0)];
        float c010 = phi_b[idxCell(i0, j0 + 1, k0)];
        float c110 = phi_b[idxCell(i0 + 1, j0 + 1, k0)];
        float c001 = phi_b[idxCell(i0, j0, k0 + 1)];
        float c101 = phi_b[idxCell(i0 + 1, j0, k0 + 1)];
        float c011 = phi_b[idxCell(i0, j0 + 1, k0 + 1)];
        float c111 = phi_b[idxCell(i0 + 1, j0 + 1, k0 + 1)];
        // tri-linear interpolate φ_b
        float c00 = lerp(c000, c100, fracx);
        float c10 = lerp(c010, c110, fracx);
        float c01 = lerp(c001, c101, fracx);
        float c11 = lerp(c011, c111, fracx);
        float c0 = lerp(c00, c10, fracy);
        float c1 = lerp(c01, c11, fracy);
        float phi_b_val = lerp(c0, c1, fracz);
        float phi_w_val = 1.0f - phi_b_val;
        if (phi_w_val < 0) phi_w_val = 0;
        return phi_w_val;
    }

    // Solve variable-density pressure Poisson to enforce incompressibility
    void projectPressure(float dt) {
        // Build RHS divergence: water divergence + bubble volume change
        std::fill(divergence.begin(), divergence.end(), 0.0f);
        // Compute water velocity divergence at cell centers
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    // outflow = (u_right - u_left + v_top - v_bottom + w_front - w_back)
                    float div = 0.0f;
                    div += u[idxU(i + 1, j, k)] - u[idxU(i, j, k)];
                    div += v[idxV(i, j + 1, k)] - v[idxV(i, j, k)];
                    div += w[idxW(i, j, k + 1)] - w[idxW(i, j, k)];
                    divergence[idxCell(i, j, k)] = div / dx; // volume flux per volume
                }
            }
        }
        // Add bubble volume change term: -(φ_b_new - φ_b_old)/dt
        for (int idx = 0; idx < Nx * Ny * Nz; ++idx) {
            float bubbleChange = (phi_b[idx] - phi_b_prev[idx]) / dt;
            divergence[idx] -= bubbleChange;
        }
        // Solve ∇·( (1/ρ) ∇ p ) = divergence/dt  (approximately)
        // We use Gauss-Seidel iterations for simplicity
        std::fill(pressure.begin(), pressure.end(), 0.0f);
        const int maxIters = 80;
        for (int iter = 0; iter < maxIters; ++iter) {
            // Sweep through cells
            for (int k = 0; k < Nz; ++k) {
                for (int j = 0; j < Ny; ++j) {
                    for (int i = 0; i < Nx; ++i) {
                        // Compute coefficient sum (diag) and neighbor contributions
                        float coefSum = 0.0f;
                        float neighborTerm = 0.0f;
                        // X neighbors
                        if (i > 0) {
                            float rho_face = effectiveRhoOnFaceX(i, j, k);
                            float coef = 1.0f / rho_face;
                            coefSum += coef;
                            neighborTerm += coef * pressure[idxCell(i - 1, j, k)];
                        }
                        if (i < Nx - 1) {
                            float rho_face = effectiveRhoOnFaceX(i + 1, j, k);
                            float coef = 1.0f / rho_face;
                            coefSum += coef;
                            neighborTerm += coef * pressure[idxCell(i + 1, j, k)];
                        }
                        // Y neighbors
                        if (j > 0) {
                            float rho_face = effectiveRhoOnFaceY(i, j, k);
                            float coef = 1.0f / rho_face;
                            coefSum += coef;
                            neighborTerm += coef * pressure[idxCell(i, j - 1, k)];
                        }
                        if (j < Ny - 1) {
                            float rho_face = effectiveRhoOnFaceY(i, j + 1, k);
                            float coef = 1.0f / rho_face;
                            coefSum += coef;
                            neighborTerm += coef * pressure[idxCell(i, j + 1, k)];
                        }
                        // Z neighbors
                        if (k > 0) {
                            float rho_face = effectiveRhoOnFaceZ(i, j, k);
                            float coef = 1.0f / rho_face;
                            coefSum += coef;
                            neighborTerm += coef * pressure[idxCell(i, j, k - 1)];
                        }
                        if (k < Nz - 1) {
                            float rho_face = effectiveRhoOnFaceZ(i, j, k + 1);
                            float coef = 1.0f / rho_face;
                            coefSum += coef;
                            neighborTerm += coef * pressure[idxCell(i, j, k + 1)];
                        }
                        if (coefSum == 0) continue; // isolated cell (should not happen unless out of domain)
                        // Right-hand side term = divergence (units 1/s)
                        float b = divergence[idxCell(i, j, k)] / dt;
                        // Gauss-Seidel update: p_new = (neighborTerm - b) / coefSum
                        float newP = (neighborTerm - b) / coefSum;
                        // Under-relaxation to improve stability (optional)
                        pressure[idxCell(i, j, k)] = newP;
                    }
                }
            }
        }
        // Correct velocities using pressure
        // interior faces
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 1; i < Nx; ++i) {
                    // X-face between cell (i-1) and (i)
                    float rho_face = effectiveRhoOnFaceX(i, j, k);
                    float pL = pressure[idxCell(i - 1, j, k)];
                    float pR = pressure[idxCell(i, j, k)];
                    // u update: Δu = -(dt/ρ_face)*(pR - pL)/dx. Here dx is cell size
                    u[idxU(i, j, k)] -= (dt / rho_face) * ((pR - pL) / dx);
                }
            }
        }
        for (int k = 0; k < Nz; ++k) {
            for (int j = 1; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    float rho_face = effectiveRhoOnFaceY(i, j, k);
                    float pB = pressure[idxCell(i, j - 1, k)];
                    float pT = pressure[idxCell(i, j, k)];
                    v[idxV(i, j, k)] -= (dt / rho_face) * ((pT - pB) / dx);
                }
            }
        }
        for (int k = 1; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    float rho_face = effectiveRhoOnFaceZ(i, j, k);
                    float pF = pressure[idxCell(i, j, k - 1)];
                    float pBa = pressure[idxCell(i, j, k)];
                    w[idxW(i, j, k)] -= (dt / rho_face) * ((pBa - pF) / dx);
                }
            }
        }
        // Enforce boundary conditions: no flow at solid boundaries (set boundary face velocities to 0)
        // X boundaries
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                u[idxU(0, j, k)] = 0.0f;
                u[idxU(Nx, j, k)] = 0.0f;
            }
        }
        // Y boundaries
        for (int i = 0; i < Nx; ++i) {
            for (int k = 0; k < Nz; ++k) {
                v[idxV(i, 0, k)] = 0.0f;
                v[idxV(i, Ny, k)] = 0.0f;
            }
        }
        // Z boundaries
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                w[idxW(i, j, 0)] = 0.0f;
                w[idxW(i, j, Nz)] = 0.0f;
            }
        }
    }

    // Advance the simulation by one time step
    void step(float dt) {
        // 1. Advect bubbles
        advectBubbles(dt);
        // 2. Update volume fractions
        updatePhiB();
        // 3. Apply forces (gravity, buoyancy, drag, surface tension)
        applyForces(dt);
        // 4. Pressure projection for incompressibility
        projectPressure(dt);
        // (Water velocity field is now divergence-free relative to mixture; bubble velocities updated from forces already)
    }

private:
    // Linear interpolation helper
    inline float lerp(float a, float b, float t) const { return a + t * (b - a); }
    // Effective density on face between cell (i-1) and i in X direction
    float effectiveRhoOnFaceX(int i, int j, int k) const {
        // face i separates cell i-1 and i
        if (i <= 0 || i >= Nx) {
            return rho_w; // outside domain, not used
        }
        int cL = idxCell(i - 1, j, k);
        int cR = idxCell(i, j, k);
        // linear blend of densities weighted by volume fractions
        float phi_b_L = phi_b[cL];
        float phi_b_R = phi_b[cR];
        float phi_b_face = 0.5f * (phi_b_L + phi_b_R);
        float phi_w_face = 1.0f - phi_b_face;
        float rho_eff = phi_b_face * rho_b + phi_w_face * rho_w;
        return rho_eff;
    }
    float effectiveRhoOnFaceY(int i, int j, int k) const {
        if (j <= 0 || j >= Ny) {
            return rho_w;
        }
        int cB = idxCell(i, j - 1, k);
        int cT = idxCell(i, j, k);
        float phi_b_B = phi_b[cB];
        float phi_b_T = phi_b[cT];
        float phi_b_face = 0.5f * (phi_b_B + phi_b_T);
        float phi_w_face = 1.0f - phi_b_face;
        float rho_eff = phi_b_face * rho_b + phi_w_face * rho_w;
        return rho_eff;
    }
    float effectiveRhoOnFaceZ(int i, int j, int k) const {
        if (k <= 0 || k >= Nz) {
            return rho_w;
        }
        int cF = idxCell(i, j, k - 1);
        int cB = idxCell(i, j, k);
        float phi_b_F = phi_b[cF];
        float phi_b_B = phi_b[cB];
        float phi_b_face = 0.5f * (phi_b_F + phi_b_B);
        float phi_w_face = 1.0f - phi_b_face;
        float rho_eff = phi_b_face * rho_b + phi_w_face * rho_w;
        return rho_eff;
    }
};

#endif // SIMULATION_H
