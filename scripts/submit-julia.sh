#!/usr/bin/env bash
#SBATCH --account=S2538
#SBATCH --time=11:59:00

set -euxo pipefail

echo "Current directory: $(pwd)"
echo "Julia binary: $(which julia)"
echo "Julia version: $(julia --version)"
echo "Script to run: $1"

julia -t 48 $1
