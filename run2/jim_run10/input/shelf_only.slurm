#!/bin/csh
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p sched_mit_hill,newnodes
#SBATCH --time=12:00:00
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dngoldberg@gmail.com
#SBATCH -o /home/dgoldber/slurm_submissions/job_%j.out
#SBATCH -e /home/dgoldber/slurm_submissions/job_%j.err
./etc/profile.d/modules.sh
module add gcc
module add mvapich2/gcc

# Initialise environment module

setenv LD_LIBRARY_PATH "$LD_LIBRARY_PATH":/home/dgoldber/petsc-3.3-p6/arch-linux2-c-opt/lib

set proj = coupled
set exp = shelf_only1km

set rootdir   = /home/dgoldber/MITgcm
set projdir   = ${rootdir}/${proj}
set outputdir = ${projdir}/${exp}/run
set inputdir  = ${projdir}/${exp}/input
set builddir  = /nobackup1/dgoldber/{$exp}/build
set rundir    = /nobackup1/dgoldber/{$exp}/run

echo $rundir
if (! -d $rundir) then
 mkdir -p $rundir
endif

echo $outputdir
if (! -d $outputdir) then
 mkdir -p $outputdir
endif

cd $rundir
ln -s $builddir/mitgcmuv .
ln -s $inputdir/* . 

# Run the program
mpirun -n 16 ./mitgcmuv

#cp STDOUT* $outputdir
#cp *.data $outputdir
#cp *.meta $outputdir
