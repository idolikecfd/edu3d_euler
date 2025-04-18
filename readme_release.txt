##############################################################
# Educationally-Designed Unstructured 3D (EDU3D) Code
#
#  --- EDU3D Euler
#
#  Instruction: type the following commands, or simply source
#               this file by typing in at the prompt:
#
#               source readme_release.txt
#
##############################################################


##############################################################
# Compile edu3d_euler
##############################################################

  echo "Compiling edu3d_euler......"

 #This is a debug version (more checks but slower).

  cd executable_debug
  make
  cd ../

 #This is an optimized version to be used in the tests below.

  cd executable_optimized
  make
  cd ../

##############################################################
# Run a truncation error analysis
##############################################################

 # Create a testcase directory

   cp -r test_source/mms_te  mms_te

 # Go in there.

   cd mms_te

 # Run the test. See the output on screen.

   ../executable_optimized/edu3d_euler

 # Return to the main directory.

   cd ..


##############################################################
# Run the cube-freestream case.
##############################################################

 # Create a testcase directory

   cp -r test_source/cube_freestream  testcase_cube_freestream

 # Go in there.

   cd testcase_cube_freestream

 # Run the test. See test_source/cube_freestream/readme.txt for details

   source readme.txt

 # Return to the main directory.

   cd ..

##############################################################
# Run the hemisphere-cylinder (half-geometry) case.
##############################################################

 # Create a testcase directory

   cp -r test_source/hc_half  testcase_hc_half

 # Go in there.

   cd testcase_hc_half

 # Run the test. See test_source/hc_half/readme.txt for details

   source readme.txt

 # Return to the main directory.

   cd ..

##############################################################
# Run the OM6 wing case.
##############################################################

 # Create a testcase directory

   cp -r test_source/om6  testcase_om6

 # Go in there.

   cd testcase_om6

 # Run the test. See test_source/om6/readme.txt for details

   source readme.txt

 # Return to the main directory.

   cd ..

##############################################################
# Now, go into each testcase directory and check the results.
# 1. Plot the boundary grid and solutions.
# 2. See fort.1000/2000 to see how the iteration converged.
##############################################################
