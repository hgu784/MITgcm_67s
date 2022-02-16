initial : create initial shape using ice only model
 
run : spin up 60year using coupled model
 jim_run3 : CTRL
 jim_run8 : change groundinig line (~/Figure/Fig.1)
 jim_run9 : change groundinig line (~/Figure/Fig.2)
 jim_run10 : change groundinig line (~/Figure/Fig.3)
 jim_run11 : change groundinig line (~/Figure/Fig.4)
 


run2 : sensitive experiments after spin up using coupled model
 jim_run3 : CTRL 1 km
 jim_run5 : CTRL 500 m
 jim_run6 : CTRL 250 m
 jim_run9 : change groundinig line (~/Figure/Fig.2) 250 m
 jim_r un10 :change groundinig line (~/Figure/Fig.3)  250 m

pig : jim's setup, 'STREAMICEthickFile' and 'bathyFile' change pig data

 subglacier : SHELFICEaddrunoff = .true., 
              SHELFICESubglFluxFile = 'subglflux.box',
 flux divergence : SHELFICEMassDynTendFile = 'massdyntend1.bin',

melt : ice only, change meltrate in data.obcs

pig_exf : exf & cal & obcs turn on

seasonal : seasonal change 



%How to make initial
 high resolution : change 'STREAMICEthickFile' and 'STREAMICEcalveMaskFile'
 change grounding line : change 'STREAMICEHBCyFile'

 cd input
 (in matlab) rdmds_init.m run
             gendata.m run



%How to build and run

 mkdir build
 cd build
 module load intel/2021.2.0
 module load impi/2021.2.0
 export LANG=en_US.UTF-8
 export LC_ALL=en_US.utf8
 ../../../tools/genmake2 -of ../../../tools/build_options/linux_amd64_ifort+mpi_ice_nas_tokyo3 -mpi -mods ../code/
 make depend
 export LANG=en_US.UTF-8
 export LC_ALL=en_US.utf8
 make -j 16

 mkdir run
 cd run
 cp ../build/mitgcmuv .
 cp ../input/* .
 pjsub job*.pbs
