#

 echo " ################################################################"
 echo "       Test case 02: Hemisphere cylinder"
 echo " ################################################################"

#############################################################
# 1. Generate grids for hemisphere-cylinder case
#############################################################

cd grids

 ######################################
 # Generate a fine grid
 ######################################


 echo " Compiling a grid generation code..."

 gfortran -O2 -o hcf_hc hcf_hc_v8p0.f90

 echo " Generating a grid..."

 ./hcf_hc > grid_generation_screenout.txt

 ######################################
 # Coarsen the fine grid.
 ######################################

 echo " Compiling a grid coarsening code..."

 gfortran -O2 -o hcf_coarsening hcf_coarsening_v3p3.f90

 echo " Coarsening the grid..."

 ./hcf_coarsening > grid_coarsening_screenout.txt

 ######################################
 # Back to the main directory
 ######################################

 cd ../

#############################################################
# 2. Select a grid, and copy it to the main directory.
#############################################################

cp grids/hc_tetra.3.ugrid hc_tet.ugrid

#############################################################
# 3. Run EDU3D-Euler
#
#    - Executable 'edu3d_euler' is in the top level.
#    - This will read input.nml in the current directory.
#    - It reads the included file 'hc_tet.mapbc' for BC info.
#############################################################

../executable_optimized/edu3d_euler

# The following output files will be generated:
#
# - hc_tet_boundary_tec.dat (tecplot boundary grid file)
# - fort.1000 (linear relaxation/preconditioner convergence)
# - fort.2000 (GCR linear convergence)

