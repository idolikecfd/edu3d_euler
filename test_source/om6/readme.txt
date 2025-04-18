#

 echo " ################################################################"
 echo "          Test case 03: ONERA M6 wing"
 echo " ################################################################"

#############################################################
# 1. Generate grids for OM6 case
#############################################################

cd grids

 ######################################
 # Generate a fine grid: level = 1
 ######################################

 echo " Compiling a grid generation code..."

 gfortran -O2 -o hcf_wing hcf_wing_v4p0.f90

 echo " Generating a grid..."

 ./hcf_wing > grid_generation_screenout.txt


 ######################################
 # Coarsen the fine grid: level = 2,3,4
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
#    Level-3 grid is a good size for testing.
#
#    Level :   nodes
#    -----------------
#        1 : 113,945
#        2 :  14,685
#        3 :   1,955
#        4 :     279
#    -----------------
#############################################################

cp grids/wing_tetra.3.ugrid  om6_tet.ugrid

#############################################################
# 3. Run EDU3D-Euler
#
#    - Executable 'edu3d_euler' is in the top level.
#    - This will read input.nml in the current directory.
#    - It reads the included file 'hc_tet.bcmap' for BC info.
#############################################################

../executable_optimized/edu3d_euler

# The following output files will be generated:
#
# - hc_tet_boundary_tec.dat (tecplot boundary grid file)
# - fort.1000 (linear relaxation/preconditioner convergence)
# - fort.2000 (GCR linear convergence)

