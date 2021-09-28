# ScintillatorMC

Project to mimic a scintilaltor based detector for proton radiography.

Compile by running source compile.sh

Prerequesite:
- Geant4
- ROOT

## Main Code 
The code is ran through a config file. An example command to run is

```shell script
./bin/Linux-g++/ScintillatorMC pCTconfig.txt $RANDOMNUMBER
```

where the variable $RANDOMNUMBER is to be replaced by a random number. The variables required in pCTconfig.txt are predefined, and they respectively mean:

- **NParticle**: The number of particle	per pencil beam	in the simulation
- **Energy**   : The energy of each particle. In the src/PrimaryGeneratorAction.cc, it can be changed between Mono and Spectrum	(for the latter this argument is not used).
- **Model**    : Which phantom model to	use (right now we are using LasVegas and SlantedEdge, but a whole library exist)
- **Angle**    : For tomographic reconstruction, the angle between the beam and	the phantom main axis
- **Thickness**: Represents the	thickness of the world which should be greater than the	phantom	(used to investigate scattering artifact)
- **ANumber**  : Atomic	mass of	the incoming particle (1 for proton, 4 for helium, etc..)
- **NPB**      : Number	of pencil beam in either direction, it will be squared to find the total number	of pencil beam.
- **sigmaY**   : Standard deviation of the pencil beam in the Y direction - [mm] (Reference 6.393 mm)
- **sigmaX**   : Standard deviation of the pencil beam in the X direction - [mm] (Reference 6.5524 mm)
- **sigma_AngX**   : Divergence of the pencil beam in the X direction - [mRad] (Reference 25 mRad)
- **sigma_AngY**   : Divergence of the pencil beam in the Y direction - [mRad] (Reference 25 mRad)


where you replace every variable by what you wishes (e.g. replace Energy by 200).
Except if you are using an XCAT phantom, where you need to add the CT line to define which phase to use (please refer to pCTconfig.txt for more details)

## Calibration

To perform a calibration scan, you need	to input a single pencil beam with the Empty phantom.

```shell script
./bin/Linux-g++/ScintillatorMC pCT_config.txt $RANDOMNUMBER
```
where the variable NPB is set to 1, which will position the single pencil beam in the middle of the phantom.

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

## Output data formatting
All results are saved in the output root file.
- The simulation parameters are saved in the header TTree. (NPB, SigmaY, SigmaZ, PB positions)
- 2D projection in each lateral plane (e.g. YZProj or YZProj_Q). The underscript Q stands for quenched which relates to light emission.
- The 3-D cumulative energy, light, LET and entries histogram in the scintillator.
- Front/Back represents single event proton radiograph in WET aggregated at the frontal or distal tracker. 
