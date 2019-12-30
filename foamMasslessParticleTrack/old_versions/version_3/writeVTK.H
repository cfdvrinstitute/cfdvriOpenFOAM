// writing vtk file for the particles.

word particleID(name(pcount));

// creating file pointer
autoPtr<OFstream> vtkFilePtr;
vtkFilePtr.reset(new OFstream(vtkDir/"particle_"+particleID+"_path.vtk"));

// writing vtk header
vtkFilePtr() << "# vtk DataFile Version 2.1" << endl << "particle_" << particleID << "_path" << endl << "ASCII" << endl << "DATASET POLYDATA" << endl;
vtkFilePtr() << nl;

// writing points header
vtkFilePtr() << "POINTS " << points.size() << " DOUBLE" << endl;

// writing point coordinates
forAll(points,pt)
{
  vtkFilePtr() << points[pt].x() << " " << points[pt].y() << " " << points[pt].z() << endl;
}
vtkFilePtr() << nl;

// writing lines header and list
vtkFilePtr() << "LINES 1 " << points.size() + 1 << endl; // it is LINES <no of lines> <size of list>
vtkFilePtr() << points.size() << endl;			 // <no of points> <point list .....>
forAll(points,pt)
{
  vtkFilePtr() << pt << endl;
}

Info << tab << "VTK file writen .. " << endl;
