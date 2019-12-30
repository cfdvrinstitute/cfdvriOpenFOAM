//////////////////////////////////////////////////
// creating coordinates based on cells on patch //
// using patch name			        //
//////////////////////////////////////////////////

word patchName;

// reading patch name
propDict.lookup("patchName") >> patchName;

// retriving patch id
label patchID(mesh.boundaryMesh().findPatchID(patchName));

// getting the cell ids
List<label> faceCells(mesh.boundary()[patchID].patch().faceCells());

// assembling the cell center coordinates as particle starting points
for(int i = 0; i < faceCells.size(); i++)
  {
    label pos(faceCells[i]);
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
