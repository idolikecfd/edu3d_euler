# Build the fine grid generator
add_executable(hcf_hc grids/hcf_hc_v8p0.f90 ../../src/debug_hooks.c)

# Build the fine grid coarsener
add_executable(hcf_coarsening grids/hcf_coarsening_v3p3.f90 ../../src/debug_hooks.c)

# Set up the test using the Python script from source directory with Python interpreter
add_test(NAME hc_half
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/hc_half.py
            --grid-exec $<TARGET_FILE:hcf_hc>
            --coarsen-exec $<TARGET_FILE:hcf_coarsening>
            --solver-exec $<TARGET_FILE:edu3d_euler>
            # Grid generator options
            --grid-input ${CMAKE_CURRENT_SOURCE_DIR}/grids/input.nml
            # Coarsening options
            --coarsen-input ${CMAKE_CURRENT_SOURCE_DIR}/grids/input_coarsen.nml
            -c 3
            # Solver options
            -i ${CMAKE_CURRENT_SOURCE_DIR}/input.nml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
