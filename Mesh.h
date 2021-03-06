#ifndef MESH_H
#define MESH_H

#include <vector>
#include <sstream>

class Mesh; // Forward declaration required for class Polygon to know
            // about class Mesh.

class Point {
public:
    Point(double _x, double _y) : x(_x), y(_y) {};
    friend std::ostream& operator<<(std::ostream& os, const Point& p);
    double x, y;
};

class Polygon {
public:
    Polygon(std::vector<int> const& _node_id, Mesh & _m) : node_id(_node_id), mesh(_m) {};
    std::vector<int> node_id;
    Mesh& mesh;
//    Point centroid(int p);
};

class Mesh {
public:
    Mesh(double xl, double xr, double yl, double yr, int nx, int ny);
    std::vector<Polygon> cell;
    std::vector<Point> node;
	std::vector<int> pointCellNeighbors(int p);
	std::vector<std::vector<double>> centroid(int p);
	
};

#endif //MESH_H
