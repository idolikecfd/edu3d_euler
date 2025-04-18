#

 echo " ################################################################"
 echo "        Test case 01: Free stream  in a cube"
 echo " ################################################################"

#############################################################
# 1. Generate grids for a cube free stream case.
#############################################################

cd grids

 ######################################
 # Generate a grid (3375 nodes)
 ######################################

 gfortran -O2 -o cube edu3d_tetgrid_cube_ptb_v6p1.f90
 echo " Generating a grid..."
 ./cube > grid_generation_screenout.txt 

 cd ../

#############################################################
# 2. Copy the grid to the main directory.
#############################################################

cp grids/tetgrid.ugrid  tetgrid.ugrid

#############################################################
# 3. Run EDU3D-Euler
#
#    - Executable 'edu3d_euler' is in the top level.
#    - This will read input.nml in the current directory.
#    - It reads the included file 'tetgrid.bcmap' for BC info.
#
#    - Exact solution is free stream (constants).
#    - Solver begins with perturbed solutions.
#
#############################################################

../executable_optimized//edu3d_euler

# The following output files will be generated:
#
# - tetgrid_tec.dat (tecplot boundary grid file)
# - fort.1000 (linear relaxation/preconditioner convergence)
# - fort.2000 (GCR linear convergence)

