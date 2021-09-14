export G4WORKDIR=$PWD
mkdir build
cd build

rm -rf *

cmake -DGeant4_DIR=/home/cacof1/Software/geant4.10.07/lib/Geant4-10.7.0 -DCMAKE_INSTALL_PREFIX=$G4WORKDIR ../
#cmake -DCMAKE_INSTALL_PREFIX=$G4WORKDIR ../

make -j8

make install

cd ../


