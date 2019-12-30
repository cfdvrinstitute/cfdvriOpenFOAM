#include "fvCFD.H"
#include "meshSearch.H"

// // function declaration and definitions
// word findPatch(List<word> cellPatchNames, label cellId)
// {
//   forAll(cellPatchNames,idx)
//     {
//       label patchId = mesh.boundaryMesh().findPatchID(cellPatchNames[idx]);
//       List<label> cells = mesh.boundary()[patchI].patch().faceCells();
//       bool found = (std::find(cells.begin(), cells.end(), cellId) != cells.end());

//       if (found)
// 	{
// 	  return cellPatchNames[idx];
// 	  break;
// 	}
//     }
// }

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
     "In this method, the patch name is read from the \"patchName\" entry in dictionary and the associated cell centers are taken as starting/ending points of particles\n"
     "The particles will be backtracked if the patch is an outlet and the particles will be advanceTracked if the patch is an inlet.\n"
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

  // reading the phi field value
  Info << nl << "Reading flux surface field value .. " ;
  surfaceScalarField phi
    (
     IOobject
     (
      "phi",
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
        #include "patchMethod.H" // here it will determine whether to use advanceTrack for backTrack
	break;
      }
    }

  label pcount(1);		// just count variable
  label totalParticleKilled(0);	// total number of particles killed due to exceding max time step

  // reading maximum time step for particle tracking
  label maxTimeStep;
  propDict.lookup("maxTimeStep") >> maxTimeStep;

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
  particleFilePtr() << "Particle_Number,Age,Distance_Traveled,killed,went_out,out_patch_name" << endl;

  // variables for counting particles
  label wallHitCount(0), wentOutCount(0), lastCellId(0);
  List<int> outletCellIds;
  List<word> outletPatchNames;

  // gathering all the outlet cells for detection
  Info << nl << "Collecting all the outlet cell indices .. ";
  forAll(mesh.boundaryMesh(), patchID)
    {
      if(sum(phi.boundaryField()[patchID]) > 0) // phi > 0 for all outflow going surfaces
  	{
	  outletCellIds.append(mesh.boundary()[patchID].patch().faceCells());
	  outletPatchNames.append(mesh.boundary()[patchID].name());
  	}
    }
  Info << "Done." << endl;

  // looping over particles list
  if (!backTrackFlag)
    {
      forAll(particlePositions, pos)
	{
	  // creating point vector field
	  point pnt = particlePositions[pos], newPnt;

	  scalar age(0), distance(0), iterCount(0);
	  List<point> points;

	  points.append(pnt);

	  Info << nl << "Tracking particle : " << pcount << endl;

	  // while (ms.findCell(point(pnt.x(),pnt.y(),pnt.z()),0,true) != -1)
	  while (mesh.findCell(point(pnt.x(),pnt.y(),pnt.z())) != -1 && iterCount < maxTimeStep)
	    {
	      // getting the cellId of nearby point
	      const label cellId = mesh.findCell(point(pnt.x(),pnt.y(),pnt.z()));

	      // last cell id
	      lastCellId = cellId;

	      // getting the velocity vector field at current cell
	      const vector velocity = U[cellId];

	      // calculating the timestep to be used
	      const scalar timeStep(0.5*std::cbrt(mesh.V()[cellId])/mag(velocity)); // 0.5*charLength/charVelocity

	      // computing displacement vector
	      const vector dst = velocity*timeStep; // forward trace

	      // computing new point position
	      newPnt = pnt + dst;

	      // computing distance traveled and time taken
	      distance += mag(newPnt - pnt);
	      age += timeStep;

	      // assigning back to new pnt
	      pnt = newPnt;

	      points.append(pnt);

	      iterCount++;
	    }

	  // checking whether the particle left through outlet
	  bool leftThroughOutlet = (std::find(outletCellIds.begin(), outletCellIds.end(), lastCellId) != outletCellIds.end());
	  word outPatchName;

	  if(iterCount >= maxTimeStep)
	    {
	      Info << tab << "Particle Killed! .. exceding maximum time step count." << endl;
	      totalParticleKilled++;
	      particleFilePtr() << pcount << "," << age << "," << distance << ",yes,-,-" << endl;
	    }
	  else if(leftThroughOutlet)
	    {
	      wentOutCount++;
	      forAll(outletPatchNames,idx)
		{
		  label patchId = mesh.boundaryMesh().findPatchID(outletPatchNames[idx]);
		  List<label> cells = mesh.boundary()[patchId].patch().faceCells();
		  bool found = (std::find(cells.begin(), cells.end(), lastCellId) != cells.end());

		  if (found)
		    {
		      outPatchName = outletPatchNames[idx];
		      break;
		    }
		}
	      Info << tab <<"particle left the domain through patch: "<< outPatchName << endl;
	      particleFilePtr() << pcount << "," << age << "," << distance << ",no,yes," << outPatchName << endl;
	    }
	  else
	    {
	      Info << tab <<"particle hit a wall and brought to rest .. " << endl;
	      wallHitCount++;
	      particleFilePtr() << pcount << "," << age << "," << distance << ",no,no,-" << endl;
	    }
	  Info << tab <<"Distance traveled : " << distance << " units." << endl;
	  Info << tab <<"Particle age : " << age << " units." << endl;

	  // particleFilePtr() << pcount << ", " << age << ", " << distance << endl;

#include "writeVTK.H"

	  pcount++;
	}

      Info << nl << "Total number of particles went through outlet = " << wentOutCount << endl;
      Info << nl << "Total number of particles hit by wall and brought to rest = " << wallHitCount << endl;
    }
  else
    {
      forAll(particlePositions, pos)
	{
	  // creating point vector field
	  point pnt = particlePositions[pos], newPnt;

	  scalar age(0), distance(0),iterCount(0);
	  List<point> points;

	  points.append(pnt);

	  Info << nl << "Tracking particle : " << pcount << endl;

	  // while (ms.findCell(point(pnt.x(),pnt.y(),pnt.z()),0,true) != -1)
	  while (mesh.findCell(point(pnt.x(),pnt.y(),pnt.z())) != -1 && iterCount < maxTimeStep)
	    {
	      // getting the cellId of nearby point
	      const label cellId = mesh.findCell(point(pnt.x(),pnt.y(),pnt.z()));

	      // getting the velocity vector field at current cell
	      const vector velocity = U[cellId];

	      // calculating the timestep to be used
	      const scalar timeStep(0.5*std::cbrt(mesh.V()[cellId])/mag(velocity)); // 0.5*charLength/charVelocity

	      // computing displacement vector
	      const vector dst = -velocity*timeStep; // back trace

	      // computing new point position
	      newPnt = pnt + dst;

	      // computing distance traveled and time taken
	      distance += mag(newPnt - pnt);
	      age += timeStep;

	      // assigning back to new pnt
	      pnt = newPnt;

	      points.append(pnt);

	      iterCount++;
	    }

	  if(iterCount >= maxTimeStep)
	    {
	      Info << tab << "Particle Killed! .. exceding maximum time step count." << endl;
	      totalParticleKilled++;
	    }
	  else
	    Info << tab <<"particle dead .. " << endl;
	  Info << tab <<"Distance traveled : " << distance << " units." << endl;
	  Info << tab <<"Particle age : " << age << " units." << endl;

	  particleFilePtr() << pcount << ", " << age << ", " << distance << endl;

#include "writeVTK.H"

	  pcount++;
	}
    }

  Info << nl << "Total number of particles killed due to exceeding max time step count = " << totalParticleKilled << endl;

  Info << nl << "particle's data : distance & age, are writen to the postProcessing/ directory." << endl;

  Info << nl << "End." << endl;
}
