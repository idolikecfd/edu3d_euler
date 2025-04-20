# Build the fine grid generator
add_executable(hcf_hc grids/hcf_hc_v8p0.f90)

# Build the fine grid coarsener
add_executable(hcf_coarsening grids/hcf_coarsening_v3p3.f90)

# echo " Generating a grid..."
# ./hcf_hc > grid_generation_screenout.txt
# echo " Coarsening the grid..."
# ./hcf_coarsening > grid_coarsening_screenout.txt
