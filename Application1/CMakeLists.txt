
find_package(SDL3 CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

add_executable(Application1 src/Application1.cpp)

target_link_libraries(Application1
    PRIVATE
        AlmondFEM
        SDL3::SDL3
        OpenGL::GL
        OpenGL::GLU)

install(TARGETS Application1 RUNTIME DESTINATION bin)
