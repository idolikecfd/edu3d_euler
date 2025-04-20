import argparse
import subprocess
import os
import sys

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run grid generation, coarsening, and edu3d solver')

    # Grid generator arguments
    parser.add_argument('--grid-exec', required=True, help='Path to the grid generator executable')
    parser.add_argument('--coarsen-exec', required=True, help='Path to the grid coarsening executable')

    # Grid generator specific arguments
    parser.add_argument('--grid-input', default='grids/input.nml', help='Path to grid generator input namelist file (default: grids/input.nml)')
    parser.add_argument('--airfoil-data-file', default='om6_wing_section_sharp.dat', help='Path to the airfoil data file (default: om6_wing_section_sharp.dat)')

    # Coarsening specific arguments
    parser.add_argument('--coarsen-input', default='grids/input_coarsen.nml', help='Path to coarsening input namelist file (default: grids/input_coarsen.nml)')

    # Solver arguments
    parser.add_argument('--solver-exec', required=True, help='Path to the edu3d_euler executable')
    parser.add_argument('-i', '--input', default='input.nml', help='Path to solver input namelist file (default: input.nml)')

    args = parser.parse_args()

    # Step 1: Run grid generator
    grid_cmd = [
        args.grid_exec,
        '-nml', args.grid_input,    # Pass the grid input file with -nml
        '--airfoil-data-file', args.airfoil_data_file, # Pass the airfoil data file
    ]

    print(f"Running grid generator: {' '.join(grid_cmd)}")
    try:
        result = subprocess.run(grid_cmd, check=True)
        if result.returncode != 0:
            print(f"Grid generation failed with exit code {result.returncode}")
            sys.exit(result.returncode)
        print("Grid generation completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Grid generation failed with exit code {e.returncode}")
        sys.exit(e.returncode)

    # Step 2: Run grid coarsening
    coarsen_cmd = [
        args.coarsen_exec,
        '-nml', args.coarsen_input,  # Pass the coarsening input file with -nml
    ]

    print(f"Running grid coarsening: {' '.join(coarsen_cmd)}")
    try:
        result = subprocess.run(coarsen_cmd, check=True)
        if result.returncode != 0:
            print(f"Grid coarsening failed with exit code {result.returncode}")
            sys.exit(result.returncode)
        print("Grid coarsening completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Grid coarsening failed with exit code {e.returncode}")
        sys.exit(e.returncode)

    # Step 3: Run solver
    solver_cmd = [
        args.solver_exec,
        '-i', args.input
    ]

    print(f"Running solver: {' '.join(solver_cmd)}")
    try:
        result = subprocess.run(solver_cmd, check=True)
        if result.returncode != 0:
            print(f"Solver failed with exit code {result.returncode}")
            sys.exit(result.returncode)
        print("Solver completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Solver failed with exit code {e.returncode}")
        sys.exit(e.returncode)

    return 0

if __name__ == "__main__":
    sys.exit(main())
