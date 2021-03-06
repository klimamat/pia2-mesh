 #include "Mesh.h"
#include <iostream>
#include <vector>

std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << "(" << p.x << "," << p.y << ")";
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

std::vector<std::vector<double>> Mesh::centroid(int p){
	std::vector<std::vector<double>> centroid={};
	
	for(int i=0; i<cell.size();i++){
		Polygon const& polygonTmp = cell[i];
			std::vector<int> points_IDS= polygonTmp.node_id;
						
			double x1=node[points_IDS[0]].x;
			double y1=node[points_IDS[0]].y;
			double x2=node[points_IDS[2]].x;
			double y2=node[points_IDS[2]].y;
			
			std::vector<double> coord_centroid={(x1+x2)/2.0,(y1+y2)/2.0};
						
			centroid.push_back(coord_centroid);
			
		}
		return centroid;	
};

