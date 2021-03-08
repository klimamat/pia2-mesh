#include "Mesh.h"
#include "Field.h"
#include "output.h"
#include <iostream>
#include <cstdlib>

int main(int iargc, char* iargv[]) {

    int nx = 20;
    int ny = 20;

    Mesh m(0.0,1.0,0.0,1.0,nx,ny);
    
    for (int i=0; i<m.node.size(); ++i) {
        std::cout << m.node[i] << "\n";
    }
    
    for (int i=0; i<m.cell.size(); ++i) {
        Polygon const& p = m.cell[i];
        std::cout << "[" << p.node_id[0];
        for (int j=1; j<p.node_id.size(); ++j) {
            std::cout << "," << p.node_id[j];
        }
        std::cout << "]\n";
    }
	
	  Field<double> W(m);
	
	  for (int i=0; i<m.cell.size(); ++i) {
		    W[i] = i;
        std::cout << W[i] << "\n";
    }
	
	 std::vector<int> pointCellNeighbors;
	 for (int i=0; i<m.node.size(); ++i) {
		   pointCellNeighbors = m.pointCellNeighbors(i);
		   std::cout << "Node no. " << i << " has neighbor cells " << pointCellNeighbors[0];
		   for (int j=1; j<pointCellNeighbors.size(); j++){
			     std::cout << ", " << pointCellNeighbors[j];
		   }
		   std::cout << "\n";
	  }
	  
	outputVTK("output.vtk",m,W);

	double plocha;
	for (int i=0; i<m.cell.size(); ++i) {
		Polygon p = m.cell[i];
		std::cout <<"Area of cell no. " << i << " is " << p.area() << "\n";
	}
	
    return 0;
	
	// number of nodes

int Mesh::nCellNodes(){
int Mesh::CellNodes(){
	for (int i=0; i<m.cell.size(); ++i) {
        Polygon const& p = m.cell[i];
        for (int j=1; j<p.node_id.size(); ++j) {
            p.CellNodes = p.node_id[j];
}
        nCellNodes = p.CellNodes.getCount(i);
}
return nCellNodes;
}
}	
}
