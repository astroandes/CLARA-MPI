#!/bin/sh
 
#PBS -q batch
#PBS -N clarita_test
#PBS -l mem=128mb
#PBS -l nodes=1:ppn=12
#PBS -M  j.e.forero.romero@gmail.com
#PBS -m abe

echo "HOLA"
module load rocks-openmpi
cd /lustre/home/ciencias/fisica/je.forero/CLARA-MPI/src/
echo "HOLA"
mpiexec -n 1 ./mine.x /lustre/home/ciencias/fisica/je.forero/CLARA-MPI/script/v0hom5tDust.input

