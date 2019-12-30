//////////////////////////////////////////////////
// creating coordinates based on cells on patch //
// using patch name			        //
//////////////////////////////////////////////////

word patchName;

// reading patch name
propDict.lookup("inletPatchName") >> patchName;

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
