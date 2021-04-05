#ifndef BC_H_
#define BC_H_
#include "Compressible.h"
#include "Field.h"

template <typename T>
class BC {
public:
	virtual void apply(Mesh const& m, Field<Compressible> & W) = 0;
};

class SlipWallBC : public BC<Compressible> {
public:
	const int location_id = 1;
	virtual void apply(Mesh const& m, Field<Compressible> & W);
};

#endif