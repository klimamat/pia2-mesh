#include "Mesh.h"
#include "Field.h"
#include "output.h"
#include "Compressible.h"
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

	for(auto &t : m.boundaryNodes){
		std::cout << t << " pos: " << m.node[t] << std::endl;
	}
	
	int edgeNum=3;
	for (int i=0; i<m.cell.size(); ++i) {
		Polygon p = m.cell[i];
		std::cout << i << " " <<edgeNum << " " << p.edgeLength(edgeNum) << "\n";
	}

	std::cout << m.edge.size() << " edges:\n";
	for(auto &e : m.edge) {
	  std::cout << e.n1 << " " << e.n2 << "\n";
	}
	
	Compressible U(1, 3, 4, 1e06);
	
	std::cout << "Preassure = " << U.p() << "\n";
	
    return 0;
}

