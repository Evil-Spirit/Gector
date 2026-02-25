# Tests for Gector geometric core

add_executable(test_gector
    test_main.cpp
    test_math.cpp
    test_nurbs.cpp
    test_brep.cpp
    test_sketch.cpp
    test_operations.cpp
)

target_link_libraries(test_gector PRIVATE gector)

target_compile_features(test_gector PRIVATE cxx_std_17)

add_test(NAME gector_tests COMMAND test_gector)
