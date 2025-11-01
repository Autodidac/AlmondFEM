#include <almond_fem/almond_fem.hpp>
#include <almond_fem/bubble_simulation.hpp>
#include <safe_io/utils.hpp>

#define SDL_MAIN_HANDLED
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <SDL3/SDL.h>
#include <SDL3/SDL_opengl.h>
#include <GL/glu.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

namespace
{
    using almond::fem::bubbles::BubbleSimulation;
    using almond::fem::bubbles::Vec3;

    void setup_opengl_state()
    {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_POINT_SMOOTH);
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        glClearColor(0.08f, 0.08f, 0.1f, 1.0f);
    }

    void draw_axes(const Vec3& origin, double scale)
    {
        glBegin(GL_LINES);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3f(origin.x, origin.y, origin.z);
        glVertex3f(origin.x + scale, origin.y, origin.z);
        glColor3f(0.0f, 1.0f, 0.0f);
        glVertex3f(origin.x, origin.y, origin.z);
        glVertex3f(origin.x, origin.y + scale, origin.z);
        glColor3f(0.0f, 0.6f, 1.0f);
        glVertex3f(origin.x, origin.y, origin.z);
        glVertex3f(origin.x, origin.y, origin.z + scale);
        glEnd();
    }

    void draw_domain_bounds(const BubbleSimulation& sim)
    {
        const auto& grid = sim.grid();
        const double w = grid.nx() * grid.cell_size();
        const double h = grid.ny() * grid.cell_size();
        const double d = grid.nz() * grid.cell_size();

        glColor4f(0.3f, 0.3f, 0.35f, 0.5f);
        glBegin(GL_LINE_LOOP);
        glVertex3f(0, 0, 0);
        glVertex3f(w, 0, 0);
        glVertex3f(w, 0, d);
        glVertex3f(0, 0, d);
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex3f(0, h, 0);
        glVertex3f(w, h, 0);
        glVertex3f(w, h, d);
        glVertex3f(0, h, d);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(0, 0, 0); glVertex3f(0, h, 0);
        glVertex3f(w, 0, 0); glVertex3f(w, h, 0);
        glVertex3f(w, 0, d); glVertex3f(w, h, d);
        glVertex3f(0, 0, d); glVertex3f(0, h, d);
        glEnd();
    }

    void draw_slice(const BubbleSimulation& sim)
    {
        const auto& g = sim.grid();
        if (!g.nx() || !g.ny() || !g.nz()) return;
        const size_t slice = g.nx() / 2;
        glDisable(GL_DEPTH_TEST);
        glBegin(GL_QUADS);
        for (size_t k = 0; k < g.nz(); ++k)
            for (size_t j = 0; j < g.ny(); ++j)
            {
                const size_t idx = g.index(slice, j, k);
                double gas = g.gas_fraction()[idx];
                if (gas <= 0.01) continue;
                float a = std::clamp((float)gas, 0.f, 1.f);
                float x0 = slice * g.cell_size();
                float x1 = (slice + 1) * g.cell_size();
                float y0 = j * g.cell_size();
                float y1 = (j + 1) * g.cell_size();
                float z0 = k * g.cell_size();
                float z1 = (k + 1) * g.cell_size();
                glColor4f(0.2f, 0.4f, 1.0f, a * 0.6f);
                glVertex3f(x0, y0, z0); glVertex3f(x1, y0, z0);
                glVertex3f(x1, y1, z0); glVertex3f(x0, y1, z0);
                glVertex3f(x0, y0, z1); glVertex3f(x1, y0, z1);
                glVertex3f(x1, y1, z1); glVertex3f(x0, y1, z1);
            }
        glEnd();
        glEnable(GL_DEPTH_TEST);
    }
}

int main()
{
    if (!SDL_Init(SDL_INIT_VIDEO))
    {
        safe_io::print("SDL_Init failed: {}", SDL_GetError());
        return 1;
    }

    int width = 1280, height = 720;
    SDL_PropertiesID props = SDL_CreateProperties();
    SDL_SetStringProperty(props, SDL_PROP_WINDOW_CREATE_TITLE_STRING, "AlmondFEM Bubble Simulation");
    SDL_SetNumberProperty(props, SDL_PROP_WINDOW_CREATE_WIDTH_NUMBER, width);
    SDL_SetNumberProperty(props, SDL_PROP_WINDOW_CREATE_HEIGHT_NUMBER, height);
    SDL_SetBooleanProperty(props, SDL_PROP_WINDOW_CREATE_OPENGL_BOOLEAN, true);
    SDL_SetBooleanProperty(props, SDL_PROP_WINDOW_CREATE_RESIZABLE_BOOLEAN, true);

    SDL_Window* window = SDL_CreateWindowWithProperties(props);
    SDL_DestroyProperties(props);
    if (!window)
    {
        safe_io::print("Window creation failed: {}", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_GLContext glctx = SDL_GL_CreateContext(window);
    if (!glctx)
    {
        safe_io::print("GL context creation failed: {}", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    SDL_GL_MakeCurrent(window, glctx);
    SDL_GL_SetSwapInterval(1);
    setup_opengl_state();

    BubbleSimulation sim(64, 72, 64, 0.45);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> off(-0.4, 0.4);
    std::uniform_real_distribution<double> rdist(0.35, 0.6);

    double W = sim.grid().nx() * sim.grid().cell_size();
    double H = sim.grid().ny() * sim.grid().cell_size();
    double D = sim.grid().nz() * sim.grid().cell_size();

    for (int i = 0; i < 24; ++i)
        sim.add_bubble(
            Vec3{ W * 0.5 + off(rng) * W * 0.2,
                 H * 0.15 + std::abs(off(rng)) * H * 0.1,
                 D * 0.5 + off(rng) * D * 0.2 },
            rdist(rng));

    bool running = true, paused = false, show_slice = false;
    float yaw = 0.9f, pitch = 0.45f, dist = 60.f;
    auto last = std::chrono::steady_clock::now();
    double acc = 0.0, dt = 0.008;

    safe_io::print("Controls: P pause, S slice, arrows rotate, +/- zoom, Esc exit.");

    while (running)
    {
        SDL_Event e;
        while (SDL_PollEvent(&e))
        {
            if (e.type == SDL_EVENT_QUIT) running = false;
            else if (e.type == SDL_EVENT_WINDOW_RESIZED)
            {
                width = e.window.data1;
                height = e.window.data2;
                glViewport(0, 0, width, height);
            }
            else if (e.type == SDL_EVENT_KEY_DOWN)
            {
                switch (e.key.scancode)
                {
                case SDL_SCANCODE_ESCAPE: running = false; break;
                case SDL_SCANCODE_P: paused = !paused; break;
                case SDL_SCANCODE_S: show_slice = !show_slice; break;
                case SDL_SCANCODE_LEFT: yaw -= 0.1f; break;
                case SDL_SCANCODE_RIGHT: yaw += 0.1f; break;
                case SDL_SCANCODE_UP: pitch = std::min(pitch + 0.05f, 1.2f); break;
                case SDL_SCANCODE_DOWN: pitch = std::max(pitch - 0.05f, -0.2f); break;
                case SDL_SCANCODE_EQUALS:
                case SDL_SCANCODE_KP_PLUS: dist = std::max(10.f, dist - 2.5f); break;
                case SDL_SCANCODE_MINUS:
                case SDL_SCANCODE_KP_MINUS: dist += 2.5f; break;
                }
            }
        }

        auto now = std::chrono::steady_clock::now();
        double frame = std::chrono::duration<double>(now - last).count();
        last = now;
        if (!paused)
        {
            acc += std::min(frame, 0.05);
            while (acc >= dt) { sim.step(dt); acc -= dt; }
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0, 0, width, height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double aspect = height > 0 ? (double)width / height : 1.0;
        gluPerspective(60.0, aspect, 0.1, 500.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        Vec3 c{ W * 0.5, H * 0.5, D * 0.5 };
        double cx = c.x + dist * std::cos(pitch) * std::sin(yaw);
        double cy = c.y + dist * std::sin(pitch);
        double cz = c.z + dist * std::cos(pitch) * std::cos(yaw);
        gluLookAt(cx, cy, cz, c.x, c.y, c.z, 0, 1, 0);

        draw_axes(c, 6.0);
        draw_domain_bounds(sim);
        if (show_slice) draw_slice(sim);

        glPointSize(6.f);
        glBegin(GL_POINTS);
        for (const auto& b : sim.bubbles())
        {
            float s = std::clamp((float)(b.radius * 0.8), 0.2f, 1.f);
            glColor3f(0.85f * s, 0.9f, 0.2f + 0.3f * s);
            glVertex3f(b.position.x, b.position.y, b.position.z);
        }
        glEnd();
        SDL_GL_SwapWindow(window);
    }

    SDL_GL_DestroyContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
