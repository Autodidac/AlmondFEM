
#include "simbubbles/simulation.hpp"
#define SDL_MAIN_HANDLED  // Use our own main, not SDL's
#include <SDL3/SDL.h>
#include <SDL3/SDL_opengl.h>
#include <GL/glu.h>

#include <iostream>

// Global simulation and rendering parameters
static FluidSim* sim = nullptr;
static bool paused = false;
static bool sliceView = false;
static float camAngleYaw = 0.0f;
static float camAnglePitch = 0.0f;
static float camDist = 100.0f; // camera distance from center
static int windowWidth = 800;
static int windowHeight = 600;

// Initialize OpenGL state (basic)
void initOpenGL() {
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glPointSize(8.0f); // base point size for bubbles
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
}

// Render the bubble particles as points (or spheres)
void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Setup camera projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = (float)windowWidth / (float)windowHeight;
    gluPerspective(60.0, aspect, 0.1, 1000.0);
    // Camera transform
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // Camera orbit around center (domain center roughly at (Nx*dx/2, Ny*dx/2, Nz*dx/2))
    float cx = sim->Nx * sim->dx * 0.5f;
    float cy = sim->Ny * sim->dx * 0.5f;
    float cz = sim->Nz * sim->dx * 0.5f;
    float camX = cx + camDist * cos(camAnglePitch) * sin(camAngleYaw);
    float camY = cy + camDist * sin(camAnglePitch);
    float camZ = cz + camDist * cos(camAnglePitch) * cos(camAngleYaw);
    gluLookAt(camX, camY, camZ, cx, cy, cz, 0, 1, 0);
    // Draw a simple coordinate axis at domain center (for orientation)
    glBegin(GL_LINES);
    glColor3f(1, 0, 0); glVertex3f(cx, cy, cz); glVertex3f(cx + 10, cy, cz);
    glColor3f(0, 1, 0); glVertex3f(cx, cy, cz); glVertex3f(cx, cy + 10, cz);
    glColor3f(0, 0, 1); glVertex3f(cx, cy, cz); glVertex3f(cx, cy, cz + 10);
    glEnd();
    // Draw particles
    glColor3f(0.9f, 0.9f, 0.2f); // yellowish color for bubbles
    glBegin(GL_POINTS);
    for (const Particle& p : sim->particles) {
        // If slice view enabled, filter out particles not near center plane (e.g., within one cell of center in X)
        if (sliceView) {
            float midX = sim->Nx * sim->dx * 0.5f;
            if (fabs(p.pos.x - midX) > sim->dx) continue;
        }
        // Set point size relative to bubble radius (rough approximation)
        // Since we can't easily change point size per-point in fixed pipeline, we rely on uniform glPointSize.
        // For a better approach, one could use billboards or actual sphere geometry per particle.
        // Here we assume similar radii for simplicity or update point size globally before drawing if needed.
        glVertex3f(p.pos.x, p.pos.y, p.pos.z);
    }
    glEnd();
}

// Main loop and event handling
int main(int argc, char** argv) {
    // Initialize SDL3 with video
    if (SDL_Init(SDL_INIT_VIDEO)) {
    std::cout << "initing sdl" << "\n";
    }else{
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    // Use OpenGL for rendering
    SDL_Window* window = SDL_CreateWindow("Bubble Simulation Demo",
        //SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        windowWidth, windowHeight,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    if (!window) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, 0); // 0 for compatibility profile
    if (!SDL_GL_CreateContext(window)) {
        std::cerr << "SDL_GL_CreateContext Error: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }



    // Initialize OpenGL state
    initOpenGL();

    // Create simulation (e.g. 50x50x50 domain, cell size 1.0)
    int Nx = 50, Ny = 50, Nz = 50;
    float cellSize = 1.0f;
    FluidSim simulation(Nx, Ny, Nz, cellSize);
    sim = &simulation;
    // Initialize bubbles: e.g., a cluster of small bubbles at bottom center
    for (int n = 0; n < 10; ++n) {
        float x = (Nx * cellSize) * 0.5f + ((rand() / (float)RAND_MAX) - 0.5f) * Nx * 0.2f;
        float y = (Ny * cellSize) * 0.1f + (rand() / (float)RAND_MAX) * Ny * 0.1f;
        float z = (Nz * cellSize) * 0.5f + ((rand() / (float)RAND_MAX) - 0.5f) * Nz * 0.2f;
        float r = 0.5f + (rand() / (float)RAND_MAX) * 0.2f; // random radius ~0.5
        sim->addBubble(Vec3(x, y, z), r);
    }
    // Print basic instructions
    std::cout << "Controls: P = pause, S = slice view toggle, Arrow keys = rotate/tilt camera\n";

    // Timing
    const float dt = 0.005f;
    uint64_t prevTicks = SDL_GetTicks();

    // Main loop
    bool running = true;
    while (running) {
        // Event handling
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                running = false;
            }
            else if (event.type == SDL_EVENT_WINDOW_RESIZED) {
                windowWidth = event.window.data1;
                windowHeight = event.window.data2;
                glViewport(0, 0, windowWidth, windowHeight);
            }
            else if (event.type == SDL_EVENT_KEY_DOWN) {
                SDL_Scancode sc = event.key.scancode;
                if (sc == SDL_SCANCODE_P) {
                    paused = !paused;
                }
                else if (sc == SDL_SCANCODE_S) {
                    sliceView = !sliceView;
                }
                else if (sc == SDL_SCANCODE_LEFT) {
                    camAngleYaw -= 0.1f;
                }
                else if (sc == SDL_SCANCODE_RIGHT) {
                    camAngleYaw += 0.1f;
                }
                else if (sc == SDL_SCANCODE_UP) {
                    camAnglePitch += 0.05f;
                    if (camAnglePitch > 1.5f) camAnglePitch = 1.5f;
                }
                else if (sc == SDL_SCANCODE_DOWN) {
                    camAnglePitch -= 0.05f;
                    if (camAnglePitch < -1.5f) camAnglePitch = -1.5f;
                }
            }
        }
        // Simulation update
        if (!paused) {
            // We use a fixed dt and possibly multiple substeps per frame if frame time is large
            // For simplicity, just one or two substeps per frame:
            sim->step(dt);
            sim->step(dt);
        }
        // Rendering
        renderScene();
        SDL_GL_SwapWindow(window);

        // Simple frame cap (~60 FPS)
        uint64_t now = SDL_GetTicks();
        float frameTime = (now - prevTicks) / 1000.0f;
        if (frameTime < (1.0f / 60.0f)) {
            SDL_Delay((uint32_t)(1000.0f * (1.0f / 60.0f - frameTime)));
        }
        prevTicks = SDL_GetTicks();
    }

    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
