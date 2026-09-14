#ifndef PEEC_PARAVIEW_H
#define PEEC_PARAVIEW_H

#include "mesh.h"
#include "dual_mesh.h"
#include "complex_matrix.h"

#include <stddef.h>

/*
 * Legacy VTK ASCII.
 *
 * Совместимо, в частности, со старыми версиями ParaView,
 * включая ParaView 5.4.
 *
 *
 * Поверхность:
 *
 *     point data:
 *         phi / V
 *         Q
 *         sigma
 *
 *
 * Реберная сетка:
 *
 *     cell data:
 *         I
 *         J
 *
 * где:
 *
 *     J_e = (I_e / S_e^e) l_e.
 */


/*
 * Реальное решение:
 *
 *     transient / lightning.
 */
int paraview_write_real_state(
    const char *directory,
    const char *prefix,
    size_t frame,
    double time,
    const Mesh *mesh,
    const DualMesh *dual,
    const double *node_voltage,
    const double *node_charge,
    const double *edge_current
);


/*
 * Комплексное решение:
 *
 *     harmonic scattering / RCS.
 *
 * Для комплексных величин сохраняются:
 *
 *     _re
 *     _im
 *     _abs.
 */
int paraview_write_complex_state(
    const char *directory,
    const char *prefix,
    const Mesh *mesh,
    const DualMesh *dual,
    const Complex *node_voltage,
    const Complex *node_charge,
    const Complex *edge_current
);


/*
 * PVD collection для временной серии.
 *
 * Сначала:
 *
 *     paraview_pvd_begin()
 *
 * на каждом сохраненном кадре:
 *
 *     paraview_pvd_add()
 *
 * в конце:
 *
 *     paraview_pvd_end().
 */
typedef struct
{
    void *file;
    int opened;

} ParaViewPvd;

void paraview_pvd_init(
    ParaViewPvd *pvd
);

int paraview_pvd_begin(
    ParaViewPvd *pvd,
    const char *filename
);

int paraview_pvd_add(
    ParaViewPvd *pvd,
    double time,
    const char *vtk_file
);

void paraview_pvd_end(
    ParaViewPvd *pvd
);

#endif
