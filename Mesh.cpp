#include "Mesh.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

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

    generateEdges();
    initLeftRight();
};

//random posuv uzlu o r
void Mesh::randomize(double r){
	double vysledek, puvodniX, puvodniY, noveX, noveY, randomR, randomPhi;
	for (int j=0; j<(node.size()); ++j){
		puvodniX = node[j].x;
		puvodniY = node[j].y;
		randomR = (double)rand()/RAND_MAX*r;
		randomPhi = (double)rand()/RAND_MAX*2*3.14159265359;
		noveX = puvodniX + randomR*cos(randomPhi);
		noveY = puvodniY + randomR*sin(randomPhi);
		node[j].x = noveX;
		node[j].y = noveY;
		std::cout << "Coordinates of node " << j << "were changed  "<< "X:" << puvodniX <<"->" << noveX <<"   Y:" << puvodniY <<"->" << noveY << "   (displacement = "<< pow(pow(puvodniX-noveX,2)+ pow(puvodniY-noveY,2),0.5)<< ")\n";
		//std::cout <<"old x = " puvodniX <<"\n";//", new x = " << noveX << "\n" <<"; old y = " puvodniY <<", new y = " << noveY << " => displacement = "<< pow(pow(puvodniX-noveX,2)+ pow(puvodniY-noveY,2),0.5)<< 
	}
	return;
}

std::vector<int> Mesh::pointCellNeighbors(int p) const {
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

double Polygon::area() const {
	double plocha, lsum, rsum;
	for (int j=0; j<node_id.size()-1; ++j) {                                   // potreuju cyklus od 0 do poctu nodu meho polygonu-1
		lsum = lsum + mesh.node[node_id[j]].x * mesh.node[node_id[j+1]].y;
		rsum = rsum + mesh.node[node_id[j+1]].x * mesh.node[node_id[j]].y;
	}
	plocha = std::abs (lsum + mesh.node[node_id[node_id.size()-1]].x * mesh.node[node_id[0]].y) - rsum - (mesh.node[node_id[0]].x * mesh.node[node_id[node_id.size()-1]].y);
	plocha = plocha*0.5;
	return plocha;
}

	// number of nodes
int Mesh::nCellNodes(){
	int CellNodes, TotNodes = 0;
	for (int i=0; i<cell.size(); ++i) {
        CellNodes = cell[i].node_id.size();
		TotNodes = TotNodes + CellNodes;	
}
return TotNodes;
};

Point Polygon::centroid() const {
			double x1=mesh.node[node_id[0]].x;
			double y1=mesh.node[node_id[0]].y;
			double x2=mesh.node[node_id[1]].x;
			double y2=mesh.node[node_id[1]].y;
			double x3=mesh.node[node_id[2]].x;
			double y3=mesh.node[node_id[2]].y;
			double x4=mesh.node[node_id[3]].x;
			double y4=mesh.node[node_id[3]].y;
			
			double barycenter1x=(x1+x2+x3)/3;
			double barycenter2x=(x3+x4+x1)/3;
			double barycenter1y=(y1+y2+y3)/3;
			double barycenter2y=(y3+y4+y1)/3;
			
			double area1=x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
			double area2=x3*(y4-y1)+x4*(y1-y3)+x1*(y3-y4);
									
	return {(barycenter1x*area1+barycenter2x*area2)/(area1+area2), (barycenter1y*area1+barycenter2y*area2)/(area1+area2)};
}

//test konvexnosti bunky (1 = je konvexni; 0 = neni konvexni)
bool Polygon::isConvex() const {
	bool isConvex;
	double u11, u12, u21, u22, v11, v12, v21, v22, u31, u32, v31, v32, w13, w23, w33; 
	//prvni index = 1 - tyka se cyklu pres vsechny uzly krome poslednich dvou
	//prvni index = 2 - tyka se predposledniho uzlu
	//prvni index = 3 - tyka se posledniho uzlu
	//druhy index - slozka prislusneho vektoru
	
//cyklus pres vsechny uzly krome poslednich dvou
			for(int j=0; j<(node_id.size()-2); ++j){
			u11 = mesh.node[node_id[j+1]].x - mesh.node[node_id[j]].x;
			u12 = mesh.node[node_id[j+1]].y - mesh.node[node_id[j]].y;
			
			v11 = mesh.node[node_id[j+2]].x - mesh.node[node_id[j]].x;
			v12 = mesh.node[node_id[j+2]].y - mesh.node[node_id[j]].y;
			
			w13 = u11 * v12 - v11 * u12;
		
			if(w13 < 0){		
							
			return 0;
				}
	   		};

//vektory u,v vedene z predposledniho uzlu   	
			u21 = mesh.node[node_id[node_id.size()-1]].x - mesh.node[node_id[node_id.size()-2]].x;
			u22 = mesh.node[node_id[node_id.size()-1]].y - mesh.node[node_id[node_id.size()-2]].y;
			
			v21 = mesh.node[node_id[0]].x - mesh.node[node_id[node_id.size()-1]].x;
			v22 = mesh.node[node_id[0]].y - mesh.node[node_id[node_id.size()-1]].y;
			
			w23 = u21 * v22 - v21 * u22;

//vektory u,v vedene z posledniho uzlu 			
			u31 = mesh.node[node_id[0]].x - mesh.node[node_id[node_id.size()-1]].x;
			u32 = mesh.node[node_id[0]].y - mesh.node[node_id[node_id.size()-1]].y;
			
			v31 = mesh.node[node_id[1]].x - mesh.node[node_id[node_id.size()-1]].x;
			v32 = mesh.node[node_id[1]].y - mesh.node[node_id[node_id.size()-1]].y;
			
			w33 = u31 * v32 - v31 * u32;
	    
		if((w13 <= 0) && (w23 <= 0) && (w33 <= 0)){
	    	isConvex = 0;
		}else{	
			isConvex = 1;
		};

return isConvex;
}

//Vygeneruje unikatni hrany v siti a hranicni body (snad) libovolne nestrukturovane site (muze obsahovat i diry)
void Mesh::generateEdges(){

 //anonymni funkce - hash ktery kazde hrane priradi unikatni cislo
	auto hashFnc = [](double _x1, double _y1, double _x2, double _y2){
		double _hash;
		_hash = (_x1 + _x2)*14554 + _x1*_x2 + (_y1 + _y2)*137 + _y1*_y2*0.013;
		return _hash;
	};
	//vygenerovani jednotlivych hran -> nektere budou dvakrat
	for(auto &plg : cell){
		Polygon const& plgTemp = plg;
		int numVertices = plg.node_id.size();
		for(int i = 0; i < numVertices-1; i++){
			//ulozeni jednotlivych hran bunky
			edge.push_back(Edge(plg.node_id[i],	plg.node_id[i+1],
						hashFnc(node[plg.node_id[i]].x, node[plg.node_id[i]].y, node[plg.node_id[i+1]].x, node[plg.node_id[i+1]].y ),*this));
		}
		//hrana ktera uzavira bunku
		edge.push_back(Edge(plg.node_id[numVertices-1], plg.node_id[0],
					hashFnc(node[plg.node_id[numVertices-1]].x, node[plg.node_id[numVertices-1]].y, node[plg.node_id[0]].x, node[plg.node_id[0]].y ),*this));
	}
	
	//seradi hrany podle hodnoty hashe
    std::sort(edge.begin(), edge.end());
	std::vector<Edge> edges_sorted = edge;
    auto last = std::unique(edge.begin(), edge.end());
	edge.erase(last,edge.end());
        
	//Vybere hrany na hranici a vyhodí jejich uzly (std::set odstrani duplicity)
	for(int i = 1; i < edges_sorted.size() - 1; i++){
		if(edges_sorted[i].hash != edges_sorted[i-1].hash && edges_sorted[i].hash != edges_sorted[i+1].hash){
			boundaryNodes.insert(edges_sorted[i].n1);
			boundaryNodes.insert(edges_sorted[i].n2);
		}
	}
	if(edges_sorted[edges_sorted.size()-1].hash != edges_sorted[edges_sorted.size()-2].hash){
			boundaryNodes.insert(edges_sorted[edges_sorted.size()-1].n1);
			boundaryNodes.insert(edges_sorted[edges_sorted.size()-1].n2);
	}
	if(edges_sorted[0].hash != edges_sorted[1].hash){
			boundaryNodes.insert(edges_sorted[0].n1);
			boundaryNodes.insert(edges_sorted[0].n2);
	}
	
}

//delka hrany bunky
double Polygon::edgeLength(int i) const {	
double edgeLength;
if(i==(node_id.size()-1)){	
	edgeLength = sqrt(pow((mesh.node[node_id[0]].x-mesh.node[node_id[i]].x),2)+pow((mesh.node[node_id[0]].y-mesh.node[node_id[i]].y),2));
}
else
{
	edgeLength = sqrt(pow((mesh.node[node_id[i+1]].x-mesh.node[node_id[i]].x),2)+pow((mesh.node[node_id[i+1]].y-mesh.node[node_id[i]].y),2));
}
return edgeLength;
}

Vector2D Edge::normal() const { return Vector2D(mesh.node[n1],mesh.node[n2]).normal(); }
Vector2D Edge::unitNormal() const { return Vector2D(mesh.node[n1],mesh.node[n2]).unitNormal(); }
Point Edge::center() const { return Point(0.5*(mesh.node[n1].x+mesh.node[n2].x),0.5*(mesh.node[n1].y+mesh.node[n2].y)); }

// pomocny 2D vektor (cislo_hrany,soused_1,soused_2)
std::vector<std::vector<int>> Mesh::edgeNeighbors() const {
	std::vector<std::vector<int>> edgeNeighbors;
	std::vector<int> pointCellNeighbors_n1;
	std::vector<int> pointCellNeighbors_n2;
	
	for(int i=0;i < edge.size();++i) {
		Edge const& e = edge[i];
		
		int n1 = e.n1;
		int n2 = e.n2;
		int Neighbor1 =-2;
		int Neighbor2 =-1;
		int t = -3;
		
		pointCellNeighbors_n1 = pointCellNeighbors(n1);
		pointCellNeighbors_n2 = pointCellNeighbors(n2);

		for (int j=0; j<pointCellNeighbors_n1.size();++j) {
		
			for (int k=0; k<pointCellNeighbors_n2.size();++k) {
			
				if( t != i && pointCellNeighbors_n1[j] == pointCellNeighbors_n2[k]){
				
					Neighbor1 = pointCellNeighbors_n1[j];
					t=i;
				}
				else if( pointCellNeighbors_n1[j] != Neighbor1 && pointCellNeighbors_n1[j] == pointCellNeighbors_n2[k]){
				
					Neighbor2 = pointCellNeighbors_n1[j];
				}
			}
		}
		if ( Neighbor2 == -1) {
			edgeNeighbors.push_back({i,Neighbor1});
		}
		else {
			edgeNeighbors.push_back({i,Neighbor1,Neighbor2});	//format 2Dvektoru: (cislo_hrany,soused_1,soused_2)
		}
	}
	return edgeNeighbors;
}

void Mesh::initLeftRight() {
	std::vector<std::vector<int>> en = edgeNeighbors();
		
	for(int i=0;i < edge.size();++i) {
		//if ( en[i].size() == 3)
		//	std::cout << en[i][1] << "," << en[i][2] << "\n";
		//else
		//	std::cout << en[i][1] << "\n";
		int cl = 0;
		int cr = 0;
		Edge & e = edge[i];
		
		Vector2D normal_vektor = edge[i].normal();
		int n1_id = e.n1;	//volani cisla bodu n1 hrany i
		int n2_id = e.n2;	//volani cisla bodu n2 hrany i
		int Neighbors_id = en[i][1];
		Point centroid = cell[Neighbors_id].centroid();
		Vector2D n1_c = Vector2D(node[n1_id],centroid);
		 
		if ( en[i].size() == 3){		//kdyz bude velikost radku = 3, tak:		 
			if(dot(n1_c,normal_vektor) > 0){
				cl = en[i][1];
				cr = en[i][2];
			}
			else {
				cr = en[i][1];
				cl = en[i][2];
			}
		}
		else {	// pro krajni hranu vzdy bunka vlevo 
			cl = en[i][1];
			cr = -1;
		}
		
		e.cl = cl;
		e.cr = cr;
	}
}

// funkce: soused vlevo
int Edge::left() const{
	return cl;	
}

// funkce: soused vpravo
int Edge::right() const{
	return cr;	
}