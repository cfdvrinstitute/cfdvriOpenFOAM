#include "fvCFD.H"
#include "meshSearch.H"

int main(int argc, char *argv[])
{

  // Define the help message for this application
  argList::addNote
    (
     "This application computes the massless particle tracks in an already-solved steady-state velocity field.\n"
     "It outputs the length of the track and time taken by the particle to reach outlet.\n"
     "\n"
     "Input arguments:\n"
     "----------------\n"
     " 1) already solved velocity field for the mesh in a latest time (or 0/) folder \n"
     " 2) method to be used (points, patch, backTrack) \n"
     " 3) particleTrackDict dictionary file in the system/ folder \n"
     "\n"
     "\"points\" Method:\n"
     "In this method, the coordinates list is read from the dictionary and used as starting points of particles\n"
     "\n"
     "\"patch\" Method:\n"
     "In this method, the patch name is read from the \"inletPatchName\" entry in dictionary and the associated cell centers are taken as starting points of particles\n"
     "\n"
     "\"backTrack\" Method:\n"
     "In this method, the outlet patch name is read from the  \"outletPatchName\" entry in dictionary and the associated cell centers are taken as ending points of particles and the streamlines are back traced.\n"
     "\n"
     "developed by - Ramkumar"
     );

  // receiving the method to be used
  argList::noParallel();
  argList::validArgs.append("method");

  // checking whether timeStep is provided
  Foam::argList args(argc, argv);
  if (!args.checkRootCase())
    {
      Foam::FatalError.exit();
    }

  Info << nl;

  // #include "setRootCase.H"
#include "createTime.H"

  Info << "Reading Mesh .. " << nl << endl;
#include "createMesh.H"
  Info << "Done .. " << endl;

  // reading method
  const word METHOD = args[1];
  label method(0);

  if(METHOD == "points")
    method = 0;
  else if(METHOD == "patch")
    method = 1;
  else if (METHOD == "backTrack")
    method = 2;
  else
    {
      Info << nl << "invalid method speicfied, check \"computePaticleTrack -help\" for list of methods." << endl;
      Foam::FatalError.exit();
    }

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

  // dictionary declaration for particleTrackDict
  dictionary propDict;

  IOobject dictProp
    (
     "particleTrackDict",
     mesh.time().system(),
     mesh,
     IOobject::MUST_READ
     );

  // Check the if the dictionary is present and follows the OF format
  if (!dictProp.typeHeaderOk<dictionary>(true))
    FatalErrorIn(args.executable()) << "Cannot open \"particleTrackDict\" dictionary file! "
				    << exit(FatalError);

  Info << nl << "Method implemented : " << METHOD << endl;

  propDict = IOdictionary(dictProp);

  // declaring the list
  List<point> particlePositions;

  // bool flag for backTrack
  bool backTrackFlag(false);

  // reading the particlePositions based on method implemented
  switch(method)
    {
    case 0:			// points
      {
        #include "pointsMethod.H"
	break;
      }
    case 1:			// patch
      {
        #include "patchMethod.H"
	break;
      }
    case 2:			// backtrack
      {
        #include "backTrackMethod.H"
	backTrackFlag = true;
    	break;
      }
    }

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

      // while (ms.findCell(point(pnt.x(),pnt.y(),pnt.z()),0,true) != -1)
      while (mesh.findCell(point(pnt.x(),pnt.y(),pnt.z())) != -1)
	{
	  // getting the cellId of nearby point
	  const label cellId = mesh.findCell(point(pnt.x(),pnt.y(),pnt.z()));

	  // getting the velocity vector field at current cell
	  const vector velocity = U[cellId];

	  // calculating the timestep to be used
	  const scalar timeStep(0.5*std::cbrt(mesh.V()[cellId])/mag(velocity)); // 0.5*charLength/charVelocity

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
