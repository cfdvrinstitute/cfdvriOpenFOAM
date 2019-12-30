///////////////////////////////////////////////////////
// back tracking the particles and from outlet patch //
///////////////////////////////////////////////////////

// reading outlet patch name
word outletPatchName;

propDict.lookup("outletPatchName") >> outletPatchName;

// retriving patch id
label patchID(mesh.boundaryMesh().findPatchID(outletPatchName));

// getting the cell ids
List<label> faceCells(mesh.boundary()[patchID].patch().faceCells());

// assembling the cell center coordinates as particle starting points
for(int i = 0; i < faceCells.size(); i++)
  {
    label pos(faceCells[i]);
    particlePositions.append(mesh.C()[pos]);
  }
