# ScintillatorMC

Project to mimic a scintilaltor based detector for proton radiography.

Compile by running source compile.sh

Prerequesite:
- Geant4
- ROOT

## Main Code 
An example command to run is

```shell script
./bin/Linux-g++/ScintillatorMC NParticle Energy Model Angle Thickness Thread ANumber NPB
```

The variables respectively mean:

- **NParticle**: The number of particle	per pencil beam	in the simulation
- **Energy**   : The energy of each particle. In the src/PrimaryGeneratorAction.cc, it can be changed between Mono and Spectrum	(for the latter this argument is not used).
- **Model**    : Which phantom model to	use (right now we are using LasVegas and SlantedEdge, but a whole library exist)
- **Angle**    : For tomographic reconstruction, the angle between the beam and	the phantom main axis
- **Thickness**: Represents the	thickness of the world which should be greater than the	phantom	(used to investigate scattering artifact)
- **Thread**   : Thread	number to append to the	output root file to distinguish	when running parallel job
- **ANumber**  : Atomic	mass of	the incoming particle (1 for proton, 4 for helium, etc..)
- **NPB**      : Number	of pencil beam in either direction, it will be squared to find the total number	of pencil beam.



where you replace every variable by what you wishes (e.g. replace Energy by 200).
Except if you are using an XCAT phantom, where you run with
```shell script
./bin/Linux-g++/ScintillatorMC NParticle Energy XCAT Angle Thickness Thread ANumber LungPhantom/choose_a_phase.root
```

## Calibration

To perform a calibration scan, you need	to input a single pencil beam with the Empty phantom.

```shell script
./bin/Linux-g++/ScintillatorMC NParticle Energy Empty Angle Thickness Thread ANumber 1
```
where the NPB=1 variable will position the single pencil beam in the middle of the phantom.


## Analysis

We have currently two code to perform analysis, the GetContrast.py and GetMTF.py codes located in the PythonCode/ folder.

To run them, the following command should be used.
```shell script
python PythonCode/GetMTF.py RootFileName FOV/mm
python PythonCode/GetContrast.py RootFileName FOV/mm
```

To note, the GetMTF.py should be used with the SlantedEdge phantom, and the GetContrast with the LasVegas phantom.

## Export
To export the data we have two choices, exporting the 2-D projection or exporting the 3-D dose map (for calibration)
Run either of the following:
```shell script
python PythonCode/Write2DHist2MatLab.py RootFileName
python PythonCode/Write3DHist2MatLab.py RootFileName
```
