#ifndef BC_H_
#define BC_H_
#include "Compressible.h"
#include "Field.h"

template <typename T>
class BC {
public:
	BC(std::vector<int> loc_in) : locations(loc_in) {};
	virtual ~BC() {};
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
	virtual ~SlipWallBC() {};
	virtual void apply(Mesh const& m, Field<Compressible> & W);
};

class ReservoirBC : public BC<Compressible> {
public:
	using BC::BC; // C++11 directive for inheritance of BC class constructors
	virtual ~ReservoirBC() {};
	virtual void apply(Mesh const& m, Field<Compressible> & W);
};
///////////////////////////////////////////////////////////////////////////////////////////////
class InletBC : public BC<Compressible> {
public:
	using BC::BC; // C++11 directive for inheritance of BC class constructors
	virtual ~InletBC() {};
	virtual void apply(Mesh const& m, Field<Compressible> & W);
};

class OutletBC : public BC<Compressible> {
public:
	using BC::BC; // C++11 directive for inheritance of BC class constructors
	virtual ~OutletBC() {};
	virtual void apply(Mesh const& m, Field<Compressible> & W);
};



#endif