
add_executable(Application1 src/Application1.cpp)

target_link_libraries(Application1 PRIVATE AlmondFEM)

install(TARGETS Application1 RUNTIME DESTINATION bin)
