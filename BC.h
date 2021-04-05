#ifndef BC_H_
#define BC_H_
#include "Compressible.h"
#include "Field.h"

template <typename T>
class BC {
public:
	BC(std::vector<int> loc_in) : locations(loc_in) {};
	std::vector<int> locations;
	bool isCorrectLocation(int loc) { 
		for (int l : locations) { 
			if (l==loc) return true;
		} 
		return false; 
	};
	virtual void apply(Mesh const& m, Field<Compressible> & W) = 0;
};

class SlipWallBC : public BC<Compressible> {
public:
	using BC::BC; // C++11 directive for inheritance of BC class constructors
	virtual void apply(Mesh const& m, Field<Compressible> & W);
};

#endif