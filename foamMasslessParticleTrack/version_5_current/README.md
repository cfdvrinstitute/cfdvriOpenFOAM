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
The folder **Make** and the C++ source files (*computeParticleTrack.C, pointsMethod.H, patchMethod.H & writeVTK.H*) has to be placed in a separate directory.
the command **wmake** executing in the directory will compile the application.

## Method of executing the application ##
The method of execution is simple. The syntax of execution is as follows

>$ computeParticleTrack method

The **particleTrackDict** file must be placed in system/ directory and the sample dictionary is given in the source folder itself.
There are 2 methods, **points** and **patch**. Each method requires corresponding input from the dictionary.

### developed by Ramkumar. ###

### Version details ###
version-1: with provided fix timestep.
version-2: with autoTimestep, inletPatchNames, points and backtrack.
version-3: with autTimestep, points and patch name, inlet or outlet will be decided by the application.
version-4: with autTimestep, points and patch name, inlet or outlet will be decided by the application & random points provision for patch method.
version-5: all features of version-4 along with the particle counts as how many are left the domain through outlet and how many are died due to hitting a wall.

### note for version 4 ###
If number of random points is less than or equal to 0 or greater than total number of cells on that patch, the random points count will
be taken to be the number of cells on that patch.

### note for version 5 ###
The columns on output file **postProcessing/particlesData.csv** are *Particle_Number,Age,Distance_Traveled,killed,went_out,out_patch_name*, these
columns are explained below.
| Column name       | explanation                                                                                                 |
|:-----------------:|:-----------------------------------------------------------------------------------------------------------:|
| Particle_Number   | each particle's identity number                                                                             |
| Age               | the time period spent by the particle inside fluid domain                                                   |
| Distance_Traveled | The distance traveled by the particle inside the fluid domain                                               |
| killed            | boolian, whether a particle get killed by exceeding maxTimeStep                                             |
| went_out          | boolian, whether a particle went out of domain through an outlet patch or get hit by a wall                 |
| out_patch_name    | the name of the patch through which the particle left the domain, "-" provided for killed and hit particles |
