#include "output.h"
#include <fstream>

void outputVTK(std::string filename, Mesh const& m, Field<double> const& u) {
	std::ofstream f;
	f.open(filename,std::ios::out);
	
	f << "# vtk DataFile Version 1.0\n";
	f << "unstructured grid\n";
	f << "ascii\n";
	f << "DATASET UNSTRUCTURED_GRID\n";
	
	f << "points " << m.node.size() << " float\n";
	
	for (int i=0; i<m.node.size(); ++i) {
		f << m.node[i] << "\n";
	}
	
	f << "cells " << m.cell.size() << " " << 5*m.cell.size() << "\n"; // TODO: upravit az bude funkce nCellNodes
	
	for (int i=0; i<m.cell.size(); ++i) {
		f << m.cell[i] << "\n";
	}
	
	f << "cell_types " << m.cell.size() << "\n";
	
	for (int i=0; i<m.cell.size(); ++i) {
		f << "9\n";
	}
	
	f << "CELL_DATA " << m.cell.size() << "\n";
 	f << "SCALARS u float\n"; 
	f << "LOOKUP_TABLE default\n";
	
    for (int i=0; i<m.cell.size(); ++i) {
		f << u[i] << "\n";
	}
	
	f.close();
}

void outputVTK(std::string filename, Mesh const& m, Field<Compressible> const& u) {
	std::ofstream f;
	f.open(filename,std::ios::out);
	
	f << "# vtk DataFile Version 1.0\n";
	f << "unstructured grid\n";
	f << "ascii\n";
	f << "DATASET UNSTRUCTURED_GRID\n";
	
	f << "points " << m.node.size() << " float\n";
	
	for (int i=0; i<m.node.size(); ++i) {
		f << m.node[i] << "\n";
	}
	
	f << "cells " << m.cell.size() << " " << 5*m.cell.size() << "\n"; // TODO: upravit az bude funkce nCellNodes
	
	for (int i=0; i<m.cell.size(); ++i) {
		f << m.cell[i] << "\n";
	}
	
	f << "cell_types " << m.cell.size() << "\n";
	
	for (int i=0; i<m.cell.size(); ++i) {
		f << "9\n";
	}
	
 	f << "SCALARS rho float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i=0; i<m.cell.size(); ++i) {
		f << u[i].rho << "\n";
	}

 	f << "SCALARS |u| float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i=0; i<m.cell.size(); ++i) {
		f << u[i].rhoU.norm()/u[i].rho << "\n";
	}		

 	f << "VECTORS u float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i=0; i<m.cell.size(); ++i) {
		f << u[i].rhoU.x/u[i].rho << ", " << u[i].rhoU.y/u[i].rho << ", 0" << "\n";
	}	

 	f << "SCALARS e float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i=0; i<m.cell.size(); ++i) {
		f << u[i].e << "\n";
	}

	f.close();
}
