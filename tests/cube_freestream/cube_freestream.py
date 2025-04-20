import argparse
import subprocess
import os
import sys

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run grid generation and edu3d solver')

    # Grid generator arguments
    parser.add_argument('--grid-exec', required=True, help='Path to the grid generator executable')
    parser.add_argument('-u', '--ugrid', default='tetgrid.ugrid', help='Specify output UGRID filename (default: tetgrid.ugrid)')
    parser.add_argument('-m', '--mapbc', default='tetgrid.mapbc', help='Specify output MAPBC filename (default: tetgrid.mapbc)')
    parser.add_argument('-tv', '--tecplot-volume', default='tetgrid_tecplot.dat', help='Specify Tecplot volume file (default: tetgrid_tecplot.dat)')
    parser.add_argument('-tb', '--tecplot-boundary', default='tetgrid_boundary_tecplot.dat', help='Specify Tecplot boundary file (default: tetgrid_boundary_tecplot.dat)')
    parser.add_argument('-vol', '--write-volume-file', action='store_true', help='Enable writing of Tecplot volume file (default: disabled)')

    # Solver arguments
    parser.add_argument('--solver-exec', required=True, help='Path to the edu3d_euler executable')
    parser.add_argument('-i', '--input', default='input.nml', help='Path to input namelist file (default: input.nml)')

    args = parser.parse_args()

    # Step 1: Run grid generator
    grid_cmd = [
        args.grid_exec,
        '-u', args.ugrid,
        '-m', args.mapbc,
        '-tb', args.tecplot_boundary
    ]

    # Add the tecplot volume file option if volume writing is enabled
    if args.write_volume_file:
        grid_cmd.extend(['-tv', args.tecplot_volume, '-vol'])

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

    # Check that grid file exists
    if not os.path.exists(args.ugrid):
        print(f"Error: Grid file {args.ugrid} was not generated")
        sys.exit(1)

    # Step 2: Run solver
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
