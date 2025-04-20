# Build the grid generator
add_executable(cube grids/edu3d_tetgrid_cube_ptb_v6p1.f90)

# Set up the test using the Python script from source directory with Python interpreter
add_test(NAME cube_freestream_combined
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cube_freestream.py
            --grid-exec $<TARGET_FILE:cube>
            -u ${CMAKE_CURRENT_BINARY_DIR}/tetgrid.ugrid
            -m ${CMAKE_CURRENT_BINARY_DIR}/tetgrid.mapbc
            -tb ${CMAKE_CURRENT_BINARY_DIR}/tetgrid_boundary_tecplot.dat
            -tv ${CMAKE_CURRENT_BINARY_DIR}/tetgrid_tecplot.dat
            --write-volume-file
            --solver-exec $<TARGET_FILE:edu3d_euler>
            -i ${CMAKE_CURRENT_SOURCE_DIR}/input.nml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

