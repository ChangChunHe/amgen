#!/bin/bash -l
#SBATCH -J amgen_repro
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --partition=cpuXeon6458
#SBATCH --exclude=c22n03,c22n04,c01n05,c01n04
#SBATCH --chdir=/share/home/scuthecc/surface_AM/amgen

# 4 个家族并行, 每个 16 线程
export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16

cd /share/home/scuthecc/surface_AM/amgen
for fam in Mg6Mn6Ni6 BaFe2O4 Cs3Cr2Br9 Mg2Cr2B2O7; do
  python -u examples/reproduce_families.py --only $fam > reproduce_$fam.log 2>&1 &
done
wait
echo "=== ALL FAMILIES DONE $(date) ==="
grep -h "PASS\|FAIL\|ERROR" reproduce_*.log
