#include "fvCFD.H"
#include "meshSearch.H"

int main(int argc, char *argv[])
{

  // Define the help message for this application
  argList::addNote
    (
     "This application computes the particle tracks in an already-solved velocity field.\n"
     "It outputs the length of the track and time taken by the particle to reach outlet.\n"
     "\n"
     "Input arguments:\n"
     "----------------\n"
     " 1) already solved velocity field for the mesh in a latest time (or 0/) folder \n"
     " 2) Timestep value to be used \n"
     " 3) particleCoordinates vector file in the system/ folder \n"
     "\n"
     "Note: This is for the purpose of Post-Processing only!, mainly created for Kapahi building ventilation project"
     "\n"
     "developed by - Ramkumar"
     );

  // receiving the time step value to be used
  argList::noParallel();
  argList::validArgs.append("timeStepValue");

  // checking whether timeStep is provided
  Foam::argList args(argc, argv);
  if (!args.checkRootCase())
    {
      Foam::FatalError.exit();
    }

  Info << nl;

  // reading timestep value

  // #include "setRootCase.H"
  #include "createTime.H"

  Info << "Reading Mesh .. " << nl << endl;
  #include "createMesh.H"
  Info << "Done .. " << endl;

  // getting timestep value
  const scalar timeStep(readScalar(IStringStream(args.args()[1])()));

  // declaring mesh search engine
  meshSearch ms(mesh);

  // total timesteps present in the case
  instantList Times = runTime.times();

  #include "checkTimeOptions.H"	// this gives the startTime of computation.

  // setting the simulation time
  runTime.setTime(Times.last(),0);

  // Create and input-output object - this holds the path to the dict and its name
  IOobject dictU
    (
     "U", // name of the file
     mesh.time().timeName(), // path to where the file is
     mesh, // reference to the mesh needed by the constructor
     IOobject::MUST_READ // indicate that reading this dictionary is compulsory
     );

  // checking whether U field is available for current timestep
  if(!exists(dictU.objectPath()) && !exists(dictU.objectPath()+".gz"))
    {
      Info << nl << "No U field for Time = " << runTime.timeName() << endl;
      Foam::FatalError.exit();
    }

  // reading the U field value
  Info << nl << "Reading velocity field value .. " ;
  volVectorField U
    (
     IOobject
     (
      "U",
      mesh.time().timeName(),
      mesh,
      IOobject::MUST_READ
      ),
     mesh
     );

  Info << "Done." << endl;

  // creating particlePositions dict file instance
  IOField<point> particlePositions
    (
     IOobject
     (
      "particlePositions",
      mesh.time().system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
      ),
     particlePositions
     );


  label pcount(1);		// just count variable

  // creating output directory
  fileName outputDir = mesh.time().path()/"postProcessing";
  mkDir(outputDir);

  // creating VTK directory
  fileName vtkDir = mesh.time().path()/"VTK";
  mkDir(vtkDir);

  // create pointer
  autoPtr<OFstream> particleFilePtr;
  particleFilePtr.reset(new OFstream(outputDir/"particlesData.csv"));

  // printing the header
  particleFilePtr() << "Particle_Number, Age, Distance_Traveled" << endl;

  // looping over particles list
  forAll(particlePositions, pos)
    {
      // creating point vector field
      point pnt = particlePositions[pos], newPnt;

      scalar age(0), distance(0);
      List<point> points;

      points.append(pnt);

      Info << nl << "Tracking particle : " << pcount << endl;

      while (ms.findCell(point(pnt.x(),pnt.y(),pnt.z()),0,true) != -1)
	{
	  // getting the cellId of nearby point
	  const label cellId = ms.findCell(point(pnt.x(),pnt.y(),pnt.z()),0,true);

	  // getting the velocity vector field at current cell
	  const vector velocity = U[cellId];

	  // computing displacement vector
	  const vector dst = velocity*timeStep;

	  // computing new point position
	  newPnt = pnt + dst;

	  // computing distance traveled and time taken
	  distance += mag(newPnt - pnt);
	  age += timeStep;

	  // assigning back to new pnt
	  pnt = newPnt;

	  points.append(pnt);
	}

      Info << tab <<"particle dead .. " << endl;
      Info << tab <<"Distance traveled : " << distance << " units." << endl;
      Info << tab <<"Particle age : " << age << " units." << endl;

      particleFilePtr() << pcount << ", " << age << ", " << distance << endl;

      #include "writeVTK.H"

      pcount++;
    }

  Info << nl << "particle's data : distance & age, are writen to the postProcessing/ directory." << endl;

  Info << nl << "End." << endl;
}
