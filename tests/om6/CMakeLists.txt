# Build the fine grid generator (level = 1)
add_executable(om6_wing grids/om6_wing_v4p0.f90)

# Build the fine grid coarsener (level = 2, 3, 4)
add_executable(om6_coarsening grids/om6_coarsening_v3p3.f90)

# echo " Generating a grid..."
# ./hcf_wing > grid_generation_screenout.txt

# echo " Coarsening the grid..."
# ./hcf_coarsening > grid_coarsening_screenout.txt
