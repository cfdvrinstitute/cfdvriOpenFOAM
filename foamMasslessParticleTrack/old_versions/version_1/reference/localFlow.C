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
  // snippet for receiving faceZoneName as input

  // Define the help message for this application
    argList::addNote
    (
        "This application computes flow rate passing through given volumetric Region.\n"
	"It outputs the volume flow rate and fluid changes per hour in that region.\n"
	"\n"
	"The volumetric region is defined by the cellZone covering that region.\n"
	"And for the flux calculation, the faceZone enveloping that region is needed.\n"
        "\n"
        "Input arguments:\n"
        "----------------\n"
        "  faceZoneName - name of the faceZone that covers the volumetric region\n"
	"  cellZoneName - name of the cellZone that covers the volumetric region\n"
	"\n"
	"developed by - Ramkumar"
    );

  // receiving the facezone name
  argList::noParallel();
  argList::validArgs.append("faceZoneName");
  argList::validArgs.append("cellZoneName");

  // checking whether facezone and cellZone are provided
  Foam::argList args(argc, argv);
  if (!args.checkRootCase())
    {
      Foam::FatalError.exit();
    }
// #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
#include "createTime.H"
#include "createMesh.H"

    // total timesteps present in the case
    instantList Times = runTime.times();

#include "checkTimeOptions.H"	// this gives the startTime of computation.

    // accessing the facezone and cellZone names
    const word faceZoneName = args[1];
    const word cellZoneName = args[2];

    // accessing zoneIDs
    label faceZoneID =  mesh.faceZones().findZoneID(faceZoneName);
    label cellZoneID =  mesh.cellZones().findZoneID(cellZoneName);

    if (faceZoneID == -1)		// confirming the existence of faceZone
      FatalErrorIn(args.executable()) << "No such faceZone exists with name :  "
				      << faceZoneName << exit(FatalError);

    if (cellZoneID == -1)		// confirming the existence of cellZone
      FatalErrorIn(args.executable()) << "No such cellZone exists with name :  "
				      << cellZoneName << exit(FatalError);

    const faceZone& faceZoneList = mesh.faceZones()[faceZoneID];
    const cellZone& cellZoneList = mesh.cellZones()[cellZoneID];

    Info << nl << "Name of the Face Zone : " << faceZoneName << endl;

    // computing total surface area of faceZone
    scalar faceZoneArea(0);

    forAll(faceZoneList, faceid)
      faceZoneArea += mesh.magSf()[faceZoneList[faceid]];

    Info <<"Total surface Area of faceZone : " << faceZoneArea << endl;

    // computing volume of cellZone
    scalar cellZoneVolume(0);
    forAll(cellZoneList, cellid)
      cellZoneVolume += mesh.V()[cellZoneList[cellid]];

    Info << nl << "Name of the Cell Zone : " << cellZoneName << endl;
    Info << "Volume of cellZone : " << cellZoneVolume << endl;

    // creating output directory and file
    // Create the output path directory
    fileName outputDir = mesh.time().path()/"postProcessing";
    // Creathe the directory
    mkDir(outputDir);

    // File pointer to direct the output to
    autoPtr<OFstream> outputFilePtr;
    // Open the file in the newly created directory
    outputFilePtr.reset(new OFstream(outputDir/"flowRate_"+faceZoneName+".dat"));

    // writing header information
    outputFilePtr() << "The volume flow information occuring in the cellZone named : "+cellZoneName+", is printed in this file." << endl;

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

	// reading the phi field value
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

	// // creating phi list in the order of face indices
	// List<scalar> phiList;
	// for(int id = 0; id < mesh.Cf().size(); id++)
	//   {
	//     phiList.append(phi.internalField()[id]);
	//   }
	// forAll(mesh.boundaryMesh(),patchI)
	//   {
	//     label patchSize = mesh.boundary()[patchI].Cf().size();
	//     for(int i = 0; i < patchSize; i++)
	//       {
	// 	phiList.append(phi.boundaryField()[patchI][i]);
	//       }
	//   }

	// creating phi list in the order of face indices
	Info << nl << "Creating Phi list .. ";
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

	// computing the total flux entering in the boundary
	scalar Qflux(0);

	forAll(faceZoneList, index)
	  {
	    // getting exact face ID
	    label faceID = faceZoneList[index];

	    // getting the phi value of corresponding face
	    scalar flux = phiList[faceID];
	    if(flux < 0)
	      Qflux -= flux;
	    else
	      Qflux += flux;
	  }

	Info << nl << "Total flow-rate in the domain: " << cellZoneName
	     << ", for the time: " << mesh.time().timeName() << ", is = " << Qflux/2.0 << endl;

	Info << nl << "Total Fluid Change Per Hour (FCPH) in the domain: " << cellZoneName
	     << ", for the time: " << mesh.time().timeName() << ", is = " << Qflux/2.0/cellZoneVolume*3600 << " per hour." << endl;

	outputFilePtr() << nl << "Total flow-rate in the domain: " << cellZoneName
	     << ", for the time: " << mesh.time().timeName() << ", is = " << Qflux/2.0 << endl;

	outputFilePtr() << nl << "Total Fluid Change Per Hour (FCPH) in the domain: " << cellZoneName
	     << ", for the time: " << mesh.time().timeName() << ", is = " << Qflux/2.0/cellZoneVolume*3600 << " per hour." << endl;
      }

    Info << nl << "Computed Output is printed to the file: " << nl << outputFilePtr().name() << endl;
    Info << nl << "End" << endl;
}
