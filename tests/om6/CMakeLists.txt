# Build the fine grid generator (level = 1)
add_executable(om6_wing grids/om6_wing_v4p0.f90 ../../src/debug_hooks.c)

# Build the fine grid coarsener (level = 2, 3, 4)
add_executable(om6_coarsening grids/om6_coarsening_v3p3.f90 ../../src/debug_hooks.c)

# Tecplot visualization project file
configure_file(view.lay ${CMAKE_CURRENT_BINARY_DIR}/view.lay COPYONLY)

# Set up the test using the Python script from source directory with Python interpreter
add_test(NAME om6
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/om6.py
            --grid-exec $<TARGET_FILE:om6_wing>
            --coarsen-exec $<TARGET_FILE:om6_coarsening>
            --solver-exec $<TARGET_FILE:edu3d_euler>
            # Grid generator options
            --grid-input ${CMAKE_CURRENT_SOURCE_DIR}/grids/input.nml
            --airfoil-data-file ${CMAKE_CURRENT_SOURCE_DIR}/grids/om6_wing_section_sharp.dat
            # Coarsening options
            --coarsen-input ${CMAKE_CURRENT_SOURCE_DIR}/grids/input_coarsen.nml
            # Solver options
            -i ${CMAKE_CURRENT_SOURCE_DIR}/input.nml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
