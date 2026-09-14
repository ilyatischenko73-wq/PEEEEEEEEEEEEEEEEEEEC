#ifndef PEEC_TOPOLOGY_H
#define PEEC_TOPOLOGY_H

#include "mesh.h"


/*
 * Строит глобальную топологию сетки и вычисляет
 * основные геометрические характеристики.
 */
int mesh_build_topology(Mesh *mesh);


#endif
