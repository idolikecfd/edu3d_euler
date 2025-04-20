cmake_minimum_required(VERSION 3.10)

# Project name
project(EDU3D_EULER LANGUAGES Fortran C)

option(EDU3D_EULER_ENABLE_DEBUG_HOOKS "Enable debug hooks that pause on crash for debugger attachment" ON)

# Find Python interpreter (used for tests)
find_package(Python REQUIRED COMPONENTS Interpreter)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 2003)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# For Debug configuration - Replicate the debug version settings
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics")

# For Release/Optimized configuration - Replicate the optimized version settings
# Note: currently om6 test does not permit the use of -ffast-math due to accuracy criteria
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native") # -ffast-math")

# List of source files
set(SOURCES
    src/edu3d_euler_input_parameter.f90
    src/edu3d_derivative_data_df5.f90
    src/edu3d_flux_functions_ddt.f90
    src/edu3d_euler_data_module.f90
    src/edu3d_euler_gradients.f90
    src/edu3d_euler_nceb_data.f90
    src/edu3d_euler_eqns.f90
    src/edu3d_euler_fluxjac_interface.f90
    src/edu3d_mms_euler.f90
    src/edu3d_euler_bc_states.f90
    src/edu3d_euler_residual.f90
    src/edu3d_euler_jacobian.f90
    src/edu3d_euler_linear_solver.f90
    src/edu3d_gcr.f90
    src/edu3d_euler_steady_solver.f90
    src/edu3d_euler_main.f90
    src/debug_hooks.c
)

# Configure debug hooks
if(EDU3D_EULER_ENABLE_DEBUG_HOOKS)
    add_definitions(-DEDU3D_EULER_ENABLE_DEBUG_HOOKS)
endif()

# Add executable
add_executable(edu3d_euler ${SOURCES})

# Enable testing
enable_testing()

# Add test cases
add_subdirectory(tests)

