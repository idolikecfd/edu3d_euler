# Build the grid generator
add_executable(tetgrid_cube grids/edu3d_tetgrid_cube_v4.f90 ../../src/debug_hooks.c)

# Set up the test using the Python script from source directory with Python interpreter
add_test(NAME tetgrid_cube
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tetgrid.py
            --grid-exec $<TARGET_FILE:tetgrid_cube>
            --solver-exec $<TARGET_FILE:edu3d_euler>
            # Solver options
            -i ${CMAKE_CURRENT_SOURCE_DIR}/input.nml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
