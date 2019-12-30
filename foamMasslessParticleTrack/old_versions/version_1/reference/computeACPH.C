/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


int main(int argc, char *argv[])
{

  // Define the help message for this application
    argList::addNote
    (
        "This application computes the Air Changes Per Hour for the given computational Domain.\n"
	"It creates a field value named ACPH for further post-processing.\n"
	"\n"
        "Input arguments:\n"
        "----------------\n"
        "  None is required except the already solved case directory"
	"\n"
	"Note: This is for the purpose of Post-Processing only!, mainly created for Kapahi building ventilation project"
	"\n"
	"developed by - Ramkumar"
    );
    #include "setRootCase.H"
    // #include "createtime.h"
    // #include "createmesh.h"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
#include "createTime.H"
#include "createMesh.H"

    // total timesteps present in the case
    instantList Times = runTime.times();

    #include "checkTimeOptions.H"	// this gives the startTime of computation.


    Info<< "Creating field ACPH\n" << endl;


    for(int i = startTime; i < runTime.times().size(); i++)
      {
	runTime.setTime(Times[i],i);

	// Create and input-output object - this holds the path to the dict and its name
	IOobject dictPhi
	  (
	   "phi", // name of the file
	   mesh.time().timeName(), // path to where the file is
	   mesh, // reference to the mesh needed by the constructor
	   IOobject::MUST_READ // indicate that reading this dictionary is compulsory
	   );

	// checking whether phi field is available for current timestep
	if(!exists(dictPhi.objectPath()) && !exists(dictPhi.objectPath()+".gz"))
	  {
	    Info << nl << "No phi for Time = " << runTime.timeName() << " -> Skipping" << endl;
	    continue;
	  }
	else
	  Info << nl << "found phi for Time = " << runTime.timeName() << endl;

	Info << nl << "Declaring ACPH field variable .. ";
	// declaring ACPH field
	volScalarField ACPH
	  (
	   IOobject
	   (
            "ACPH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	    ),
	   mesh,
	   dimensionedScalar("p", dimensionSet(0,0,-1,0,0,0,0), 0.) // p is for some reference bc value, just dummy
	   );

	Info << "Done." << endl;

	// reading the phi field value
	Info << nl << "reading phi value .. ";
	surfaceScalarField phi
	  (
	   IOobject
	   (
	    "phi",
	    mesh.time().timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT
	    ),
	   mesh
	   );
	Info << "Done." << endl;

	Info << nl << "Creating Phi list .. ";

	// creating phi list in the order of face indices
	List<scalar> phiList(mesh.faces().size());
	label count(0);
	for(int id = 0; id < mesh.Cf().size(); id++)
	  {
	    phiList[count] = phi.internalField()[id];
	    count++;
	  }
	forAll(mesh.boundaryMesh(),patchI)
	  {
	    label patchSize = mesh.boundary()[patchI].Cf().size();
	    for(int i = 0; i < patchSize; i++)
	      {
		phiList[count] = phi.boundaryField()[patchI][i];
		count++;
	      }
	  }
	Info << "Done .. " << endl;

	Info << nl << "Computing ACPH .. ";
	// calculating ACPH per cell
	for(label cell = 0; cell < mesh.cells().size(); cell++)
	  {
	    // computing current cell volume
	    scalar vol = mesh.V()[cell];

	    // creating faceList
	    const labelUList faceList = mesh.cells()[cell];

	    // looping through faces for computing flux
	    scalar Qflux(0);
	    forAll(faceList,face)
	      {
		// getting exact face ID
		label faceID = faceList[face];

		// getting the phi value of corresponding face
		scalar flux = phiList[faceID];
		if(flux < 0)
		  Qflux -= flux;
		else
		  Qflux += flux;
	      }

	    // computing acph
	    // scalar acph = Qflux/2.0/vol*3600.0;
	    scalar acph = Qflux/2.0/vol;

	    ACPH[cell] = acph;
	  }
	ACPH.write();

	Info << "Done .. "<< endl;
      }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
