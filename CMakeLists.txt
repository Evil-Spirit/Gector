cmake_minimum_required(VERSION 3.16)
project(Gector VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------------------------
# Compiler warnings
# ------------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
    add_compile_options(/W4)
endif()

# ------------------------------------------------------------------------------
# Library: gector
# ------------------------------------------------------------------------------
add_library(gector
    src/nurbs/nurbs_curve.cpp
    src/nurbs/nurbs_surface.cpp
    src/brep/topology.cpp
    src/brep/brep_builder.cpp
    src/sketch/sketch_entity.cpp
    src/sketch/line.cpp
    src/sketch/arc.cpp
    src/sketch/circle.cpp
    src/sketch/spline_curve.cpp
    src/sketch/sketch.cpp
    src/operations/extrusion.cpp
    src/operations/revolve.cpp
    src/operations/boolean_op.cpp
    src/operations/fillet.cpp
)

target_include_directories(gector
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------
enable_testing()
add_subdirectory(tests)

# ------------------------------------------------------------------------------
# Install rules
# ------------------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS gector
    EXPORT  GectorTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
