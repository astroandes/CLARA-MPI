#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -pe mvapich48 96
#$ -l h_rt=1:00:00
#$ -binding linear:48
#$ -m be
#$ -M j.e.forero.romero@gmail.com
#$ -q standard.q
#$ -N T5_F

source /etc/profile.d/modules.sh
module load sge
module load hydra mvapich2/gcc/64/1.6-qlc

cd /home/extforer/CLARA-MPI/src/
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v0hom5tDust.input
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v0hom5tNoDust.input
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v100hom5tDust.input
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v100hom5tNoDust.input
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v200hom5tDust.input
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v200hom5tNoDust.input
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v300hom5tDust.input
mpistart -np 96 ./mine.x /home/extforer/CLARA-MPI/script/v300hom5tNoDust.input
