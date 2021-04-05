#ifndef INIT_H
#define INIT_H
#include <memory>
#include "Mesh.h"
#include "Field.h"
#include "Compressible.h"

void initSod(Mesh *& m, Field<Compressible> *& W); 

#endif