##How to make initial(in local machine)##

 cd initial/
 mkdir build
 cd build/
 ../../../../tools/genmake2 -mods ../code -optfile ../../../../tools/build_options/linux_amd64_gfortran
 make depend
 make -j 16

 mkdir ../run
 cd ../run/
 cp ../build/mitgcmuv .
 cp ../input/* .
 ./mitgcmuv




##How to run coupled model(in remote machine)##
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 (sNx,sNy) = (30,20),
 (nPx,nPy) = (2,5) --> mpi proc=10, 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 cd coupled/input
 
 (in matlab) rdmds_init.m  run
             gendata.m  run

 mkdir ../build
 cd ../build/
 module load intel/2021.2.0
 module load impi/2021.2.0
 export LANG=en_US.UTF-8
 export LC_ALL=en_US.utf8
 ../../../../tools/genmake2 -of ../../../../tools/build_options/linux_amd64_ifort+mpi_ice_nas_tokyo3 -mpi -mods ../code/
 make depend
 export LANG=en_US.UTF-8
 export LC_ALL=en_US.utf8
 make -j 16

 mkdir ../run
 cd ../run/
 cp ../build/mitgcmuv .
 cp ../input/* .
 pjsub job*.pbs





##How to run iceshelf model(in local machine)##

 cd iceshelf/

 %  please download INPUT FILE from Google Drive.
 
 mkdir build
 cd build/
 ../../../../tools/genmake2 -mods ../code -optfile ../../../../tools/build_options/linux_amd64_gfortran
 make depend
 make -j 16

 mkdir ../run
 cd ../run/
 cp ../build/mitgcmuv .
 cp ../input/* .
 ./mitgcmuv

