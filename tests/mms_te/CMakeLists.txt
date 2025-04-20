# Build the grid generator
add_executable(tetgrid_cube grids/edu3d_tetgrid_cube_v4.f90)

# ./tetgrid_cube
# 8,8,8
# mv tetgrid.ugrid tetgrid_08x08x08.ugrid

# ./tetgrid_cube
# 16,16,16
# mv tetgrid.ugrid tetgrid_16x16x16.ugrid

# ./tetgrid_cube
# 32,32,32
# mv tetgrid.ugrid tetgrid_32x32x32.ugrid
