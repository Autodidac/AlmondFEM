#include <almond_fem/almond_fem.hpp>
#include <almond_fem/bubble_simulation.hpp>

#include <safe_io/utils.hpp>

#define SDL_MAIN_HANDLED
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
        glVertex3f(static_cast<GLfloat>(origin.x), static_cast<GLfloat>(origin.y), static_cast<GLfloat>(origin.z));
        glVertex3f(static_cast<GLfloat>(origin.x + scale), static_cast<GLfloat>(origin.y), static_cast<GLfloat>(origin.z));

        glColor3f(0.0f, 1.0f, 0.0f);
        glVertex3f(static_cast<GLfloat>(origin.x), static_cast<GLfloat>(origin.y), static_cast<GLfloat>(origin.z));
        glVertex3f(static_cast<GLfloat>(origin.x), static_cast<GLfloat>(origin.y + scale), static_cast<GLfloat>(origin.z));

        glColor3f(0.0f, 0.6f, 1.0f);
        glVertex3f(static_cast<GLfloat>(origin.x), static_cast<GLfloat>(origin.y), static_cast<GLfloat>(origin.z));
        glVertex3f(static_cast<GLfloat>(origin.x), static_cast<GLfloat>(origin.y), static_cast<GLfloat>(origin.z + scale));
        glEnd();
    }

    void draw_domain_bounds(const BubbleSimulation& simulation)
    {
        const auto& grid = simulation.grid();
        const double width = static_cast<double>(grid.nx()) * grid.cell_size();
        const double height = static_cast<double>(grid.ny()) * grid.cell_size();
        const double depth = static_cast<double>(grid.nz()) * grid.cell_size();

        glColor4f(0.3f, 0.3f, 0.35f, 0.5f);
        glBegin(GL_LINE_LOOP);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(static_cast<GLfloat>(width), 0.0f, 0.0f);
        glVertex3f(static_cast<GLfloat>(width), 0.0f, static_cast<GLfloat>(depth));
        glVertex3f(0.0f, 0.0f, static_cast<GLfloat>(depth));
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex3f(0.0f, static_cast<GLfloat>(height), 0.0f);
        glVertex3f(static_cast<GLfloat>(width), static_cast<GLfloat>(height), 0.0f);
        glVertex3f(static_cast<GLfloat>(width), static_cast<GLfloat>(height), static_cast<GLfloat>(depth));
        glVertex3f(0.0f, static_cast<GLfloat>(height), static_cast<GLfloat>(depth));
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, static_cast<GLfloat>(height), 0.0f);

        glVertex3f(static_cast<GLfloat>(width), 0.0f, 0.0f);
        glVertex3f(static_cast<GLfloat>(width), static_cast<GLfloat>(height), 0.0f);

        glVertex3f(static_cast<GLfloat>(width), 0.0f, static_cast<GLfloat>(depth));
        glVertex3f(static_cast<GLfloat>(width), static_cast<GLfloat>(height), static_cast<GLfloat>(depth));

        glVertex3f(0.0f, 0.0f, static_cast<GLfloat>(depth));
        glVertex3f(0.0f, static_cast<GLfloat>(height), static_cast<GLfloat>(depth));
        glEnd();
    }

    void draw_slice(const BubbleSimulation& simulation)
    {
        const auto& grid = simulation.grid();
        if (grid.nx() == 0 || grid.ny() == 0 || grid.nz() == 0)
        {
            return;
        }

        const std::size_t slice = grid.nx() / 2;
        glDisable(GL_DEPTH_TEST);
        glBegin(GL_QUADS);
        for (std::size_t k = 0; k < grid.nz(); ++k)
        {
            for (std::size_t j = 0; j < grid.ny(); ++j)
            {
                const std::size_t idx = grid.index(slice, j, k);
                const double gas = grid.gas_fraction()[idx];
                if (gas <= 0.01)
                {
                    continue;
                }
                const float alpha = static_cast<float>(std::clamp(gas, 0.0, 1.0));
                const float x0 = static_cast<float>(slice * grid.cell_size());
                const float x1 = static_cast<float>((slice + 1) * grid.cell_size());
                const float y0 = static_cast<float>(j * grid.cell_size());
                const float y1 = static_cast<float>((j + 1) * grid.cell_size());
                const float z0 = static_cast<float>(k * grid.cell_size());
                const float z1 = static_cast<float>((k + 1) * grid.cell_size());
                glColor4f(0.2f, 0.4f, 1.0f, alpha * 0.6f);
                glVertex3f(x0, y0, z0);
                glVertex3f(x1, y0, z0);
                glVertex3f(x1, y1, z0);
                glVertex3f(x0, y1, z0);

                glVertex3f(x0, y0, z1);
                glVertex3f(x1, y0, z1);
                glVertex3f(x1, y1, z1);
                glVertex3f(x0, y1, z1);
            }
        }
        glEnd();
        glEnable(GL_DEPTH_TEST);
    }
}

int main()
{
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        safe_io::print("SDL_Init failed: {}", SDL_GetError());
        return 1;
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

    int width = 1280;
    int height = 720;

    SDL_Window* window = SDL_CreateWindow("AlmondFEM Bubble Simulation",
        width,
        height,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    if (!window)
    {
        safe_io::print("SDL_CreateWindow failed: {}", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_GLContext context = SDL_GL_CreateContext(window);
    if (!context)
    {
        safe_io::print("SDL_GL_CreateContext failed: {}", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    SDL_GL_MakeCurrent(window, context);
    SDL_GL_SetSwapInterval(1);

    setup_opengl_state();

    BubbleSimulation simulation(64, 72, 64, 0.45);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> offset(-0.4, 0.4);
    std::uniform_real_distribution<double> radius_dist(0.35, 0.6);

    const double domain_width = static_cast<double>(simulation.grid().nx()) * simulation.grid().cell_size();
    const double domain_height = static_cast<double>(simulation.grid().ny()) * simulation.grid().cell_size();
    const double domain_depth = static_cast<double>(simulation.grid().nz()) * simulation.grid().cell_size();

    for (int i = 0; i < 24; ++i)
    {
        const double base_x = domain_width * 0.5 + offset(rng) * domain_width * 0.2;
        const double base_y = domain_height * 0.15 + std::abs(offset(rng)) * domain_height * 0.1;
        const double base_z = domain_depth * 0.5 + offset(rng) * domain_depth * 0.2;
        const double radius = radius_dist(rng);
        simulation.add_bubble(Vec3{base_x, base_y, base_z}, radius);
    }

    bool running = true;
    bool paused = false;
    bool show_slice = false;

    float yaw = 0.9f;
    float pitch = 0.45f;
    float distance = 60.0f;

    auto last_time = std::chrono::steady_clock::now();
    double accumulator = 0.0;
    const double fixed_dt = 0.008;

    safe_io::print("Controls: P pause/resume, S toggle volume slice, arrow keys rotate camera, +/- adjust zoom, Esc exit.");

    while (running)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_EVENT_QUIT)
            {
                running = false;
            }
            else if (event.type == SDL_EVENT_WINDOW_RESIZED)
            {
                width = event.window.data1;
                height = event.window.data2;
                glViewport(0, 0, width, height);
            }
            else if (event.type == SDL_EVENT_KEY_DOWN)
            {
                switch (event.key.keysym.sym)
                {
                case SDLK_ESCAPE:
                    running = false;
                    break;
                case SDLK_p:
                    paused = !paused;
                    break;
                case SDLK_s:
                    show_slice = !show_slice;
                    break;
                case SDLK_LEFT:
                    yaw -= 0.1f;
                    break;
                case SDLK_RIGHT:
                    yaw += 0.1f;
                    break;
                case SDLK_UP:
                    pitch = std::min(pitch + 0.05f, 1.2f);
                    break;
                case SDLK_DOWN:
                    pitch = std::max(pitch - 0.05f, -0.2f);
                    break;
                case SDLK_EQUALS:
                case SDLK_PLUS:
                    distance = std::max(10.0f, distance - 2.5f);
                    break;
                case SDLK_MINUS:
                case SDLK_UNDERSCORE:
                    distance += 2.5f;
                    break;
                default:
                    break;
                }
            }
        }

        const auto now = std::chrono::steady_clock::now();
        const double frame_dt = std::chrono::duration<double>(now - last_time).count();
        last_time = now;

        if (!paused)
        {
            accumulator += std::min(frame_dt, 0.05);
            while (accumulator >= fixed_dt)
            {
                simulation.step(fixed_dt);
                accumulator -= fixed_dt;
            }
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0, 0, width, height);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        const double aspect = height > 0 ? static_cast<double>(width) / static_cast<double>(height) : 1.0;
        gluPerspective(60.0, aspect, 0.1, 500.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        const Vec3 center{
            domain_width * 0.5,
            domain_height * 0.5,
            domain_depth * 0.5,
        };

        const double cam_x = center.x + distance * std::cos(pitch) * std::sin(yaw);
        const double cam_y = center.y + distance * std::sin(pitch);
        const double cam_z = center.z + distance * std::cos(pitch) * std::cos(yaw);
        gluLookAt(cam_x, cam_y, cam_z, center.x, center.y, center.z, 0.0, 1.0, 0.0);

        draw_axes(center, 6.0);
        draw_domain_bounds(simulation);
        if (show_slice)
        {
            draw_slice(simulation);
        }

        glPointSize(6.0f);
        glBegin(GL_POINTS);
        for (const auto& bubble : simulation.bubbles())
        {
            const float color_scale = static_cast<float>(std::clamp(bubble.radius * 0.8, 0.2, 1.0));
            glColor3f(0.85f * color_scale, 0.9f, 0.2f + 0.3f * color_scale);
            glVertex3f(static_cast<GLfloat>(bubble.position.x), static_cast<GLfloat>(bubble.position.y), static_cast<GLfloat>(bubble.position.z));
        }
        glEnd();

        SDL_GL_SwapWindow(window);
    }

    SDL_GL_DestroyContext(context);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
