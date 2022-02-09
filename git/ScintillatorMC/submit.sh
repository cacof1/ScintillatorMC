#!/bin/bash -l
#$ -l h_rt=12:00:0
#$ -l mem=5G 
#S -wd /home/rmapcco/Software/ProtonMC 
source /home/rmapcco/Software/root/bin/thisroot.sh 
source /home/rmapcco/Software/geant4_10_01_p03/share/Geant4-10.1.3/geant4make/geant4make.sh
cd /home/rmapcco/Software/ScintillatorMC
hadd LasVegas_200.root LasVegas_200_0.0_*.root


