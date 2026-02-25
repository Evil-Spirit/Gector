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

# ---------------------------------------------------------------------------
# OBJ visual-export tests
# ---------------------------------------------------------------------------
add_executable(test_obj_export
    test_main.cpp
    test_obj_export.cpp
)

target_link_libraries(test_obj_export PRIVATE gector)
target_compile_features(test_obj_export PRIVATE cxx_std_17)

# Create output directory alongside the binary at build time
add_custom_command(TARGET test_obj_export POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_FILE_DIR:test_obj_export>/obj_output"
    COMMENT "Creating obj_output directory for visual test output"
)

add_test(NAME obj_export_tests
    COMMAND test_obj_export
    WORKING_DIRECTORY $<TARGET_FILE_DIR:test_obj_export>
)
