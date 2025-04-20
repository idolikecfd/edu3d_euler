import argparse
import subprocess
import os
import sys

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run grid generation and edu3d solver')

    # Grid generator arguments
    parser.add_argument('--grid-exec', required=True, help='Path to the grid generator executable')

    # Solver arguments
    parser.add_argument('--solver-exec', required=True, help='Path to the edu3d_euler executable')
    parser.add_argument('-i', '--input', default='input.nml', help='Path to solver input namelist file (default: input.nml)')

    args = parser.parse_args()

    # Step 1: Run 8x8x8 grid generator
    grid_cmd = [
        args.grid_exec,
        '--x-size', '8',
        '--y-size', '8',
        '--z-size', '8',
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

    # Step 2: Run 16x16x16 grid generator
    grid_cmd = [
        args.grid_exec,
        '--x-size', '16',
        '--y-size', '16',
        '--z-size', '16',
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

    # Step 3: Run 32x32x32 grid generator
    grid_cmd = [
        args.grid_exec,
        '--x-size', '32',
        '--y-size', '32',
        '--z-size', '32',
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

    # Step 4: Run solver
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
