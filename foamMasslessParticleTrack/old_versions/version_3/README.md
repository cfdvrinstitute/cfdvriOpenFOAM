# Custom Application for tracking the massless particles in OpenFOAM #
the application **computeParticleTrack** is built to trace the massless particles that can be made to go along the *already solved* flow field.
the main motive of creating this appliction is to compute the residence time and distance traveled by the particles.

## Existing computing mechanisms in OpenFoam ##
There are few existing particle tracking applications in OpenFOAM, but none of them are massless. For example, the **icoUncoupledKinematicParcelFoam** is a particle tracking application but it involves the real particle tracking stuff like mass, drag, gravity and collisions.
Hence the custom application is made.

## Details on custom application ##
Name: computeParticleTrack
This application solves for massless particle tracks on already solved steady state flow field. The particle coordinates, inlet patch name and outlet patch name can be given in the dictionary file called **particleTrackDict**.
The outputs from the applicaiton are as follows,
  * the residence time of particle (*age* in this application).
  * the distance traveled by the particle in the domain.
  * VTK file for each particle that will show the path of travel by the particle in the domain.

The residence time and distance travelled are printed the terminal in addition to writing in the file **particlesData.csv** under **postProcessing/** folder. The vtk file for each particle is created in the **VTK/** folder in the case directory.

## compiling application ##
The folder **Make** and the C++ files (*computeParticleTrack.C, pointsMethod.H, patchMethod.H, backTrackMethod.H & writeVTK.H*) has to be placed in a separate directory.
the command **wmake** executing in the directory will compile the application.

## Method of executing the application ##
The method of execution is simple. The syntax of execution is as follows

>$ computeParticleTrack method

The **particleTrackDict** file must be placed in system/ directory and the sample dictionary is given in the source folder itself.
There are 3 methods, **points**, **patch** and **backTrack**. Each method requires corresponding input from the dictionary.

### developed by Ramkumar. ###

### Version details ###
version-1: with provided fix timestep.
version-2: with autoTimestep, inletPatchNames, points and backtrack.
version-3: with autTimestep, points and patch name, inlet or outlet will be decided by the application.
