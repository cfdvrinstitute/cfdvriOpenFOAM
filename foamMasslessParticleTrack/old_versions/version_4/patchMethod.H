//////////////////////////////////////////////////
// creating coordinates based on cells on patch //
// using patch name			        //
//////////////////////////////////////////////////


// reading patch name
word patchName;
propDict.lookup("patchName") >> patchName;

// reading number of random points to be taken on patch
label nRandomPoints;
propDict.lookup("nRandomPoints") >> nRandomPoints;

// retriving patch id
label patchID(mesh.boundaryMesh().findPatchID(patchName));

// getting the cell ids
List<label> faceCells(mesh.boundary()[patchID].patch().faceCells());

// declaring chosen cells list
List<label> chosenNums, chosenCells;

label count(0), randNum, cellCount(faceCells.size());

bool foundFlag;

if(nRandomPoints <= 0 || nRandomPoints >= cellCount)
  chosenCells = faceCells;
 else
   {
     while(count < nRandomPoints)
       {
	 randNum = std::rand()%(cellCount);

	 foundFlag = false;

	 for(int i = 0; i < chosenNums.size(); i++)
	   {
	     if(randNum == chosenNums[i])
	       {
		 foundFlag = true;
		 break;
	       }
	   }

	 if(foundFlag)
	   continue;
	 else
	   {
	     chosenCells.append(faceCells[randNum]);
	     count++;
	   }
       }
   }

// assembling the cell center coordinates as particle starting points
for(int i = 0; i < chosenCells.size(); i++)
  {
    label pos(chosenCells[i]);
    particlePositions.append(mesh.C()[pos]);
  }

// reading phi field value
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

// determining inlet or outlet patch
if(sum(phi.boundaryField()[patchID]) < 0)
  Info << nl << "Patch: " << patchName << " is of InFlow type.. Enabling advanceTrack.." << endl;
 else
   {
     Info << nl << "Patch: " << patchName << " is of OutFlow type.. Enabling backTrack.." << endl;
     backTrackFlag = true;
   }
