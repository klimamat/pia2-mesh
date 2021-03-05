#include "Mesh.h"
#include <iostream>
#include <vector>
#include <cmath>
std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << p.x << " " << p.y << " 0.0";
    return os;
};

// vypis geometrie polygonu
std::ostream& operator<<(std::ostream& os, const Polygon& p) {
	os << p.node_id.size();
	for (int j=p.node_id.size()-1; j>=0; j--) {
            os << " " << p.node_id[j];
        }
    return os;
};

Mesh::Mesh(double xl, double xr, double yl, double yr, int nx, int ny) {
    double dx = (xr - xl)/double(nx);
    double dy = (yr - yl)/double(ny);
    
    // Create nodes
    for (int i=0; i<nx+1; ++i) {
        for (int j=0; j<ny+1; ++j) {
            node.push_back({xl + i*dx, yl + j*dy});
        }
    }
    
    // Create cells
    for (int i=0; i<nx; ++i) {
        for (int j=0; j<ny; ++j) {
            cell.push_back(Polygon({i*(ny+1) + j+1,
                                    i*(ny+1) + j,
                                    (i+1)*(ny+1) + j,
                                    (i+1)*(ny+1) + j + 1},*this));
        }
    }
};

std::vector<int> Mesh::pointCellNeighbors(int p){
	std::vector<int> pointCellNeighbors;
	
	if (p < 0 || p >= node.size()){
		std::cout << "Warning: selected point does not exist.";
	}
	
	for(int i=0; i<cell.size();i++){
		Polygon const& polygonTmp = cell[i];
		for(int j=0; j<polygonTmp.node_id.size(); j++){
			if (polygonTmp.node_id[j]==p){
				pointCellNeighbors.push_back(i);
			}
		}
	}
	return pointCellNeighbors;
};

double Polygon::area(){
	double plocha, lsum, rsum;
	for (int j=0; j<(node_id.size()-1.0); ++j) {                                   // potreuju cyklus od 0 do poctu nodu meho polygonu-1
		lsum = lsum + mesh.node[node_id[j]].x * mesh.node[node_id[j+1]].y;
		rsum = rsum + mesh.node[node_id[j+1]].x * mesh.node[node_id[j]].y;
	}
	plocha = std::abs (lsum + mesh.node[node_id[node_id.size()-1.0]].x * mesh.node[node_id[0]].y) - rsum - (mesh.node[node_id[0]].x * mesh.node[node_id[node_id.size()-1.0]].y);
	plocha = plocha*0.5;
	return plocha;
}

// number of nodes

int Mesh::nCellNodes();
	for (int i=0; i<node.size();i++){
 	nCellNodes = ncellNodes.getCount(i);
  	cout << "Sum of nodes is" << nCellNodes();
    
return nCellNodes;
