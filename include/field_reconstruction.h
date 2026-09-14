#ifndef PEEC_FIELD_RECONSTRUCTION_H
#define PEEC_FIELD_RECONSTRUCTION_H

#include "mesh.h"
#include "dual_mesh.h"

#include <stddef.h>


/*
 * Восстановление узловых зарядов:
 *
 *     P Q = V.
 *
 * LU(P) желательно подготовить заранее.
 */
int field_reconstruct_charge(
    size_t n_nodes,
    const double *P,
    const double *node_voltage,
    double *node_charge
);


/*
 * Поверхностная плотность заряда:
 *
 *     sigma_j = Q_j / S_j^v.
 */
int field_compute_surface_charge_density(
    const DualMesh *dual_mesh,
    const double *node_charge,
    double *surface_charge_density
);


/*
 * Поверхностный ток:
 *
 *     J_e =
 *
 *     (I_e / w_e) e_e.
 *
 * Размер массива:
 *
 *     3 * Ne.
 */
int field_compute_surface_current_density(
    const Mesh *mesh,
    const DualMesh *dual_mesh,
    const double *edge_current,
    double *surface_current_density
);


/*
 * Скалярный потенциал в произвольной точке пространства.
 */
int field_compute_potential_at_point(
    const DualMesh *dual_mesh,
    const double *node_charge,
    const double point[3],
    int quadrature_order,
    double *phi
);


/*
 * Электрическое поле в произвольной точке.
 */
int field_compute_electric_field_at_point(
    const Mesh *mesh,
    const DualMesh *dual_mesh,
    const double *node_charge,
    const double *edge_current_derivative,
    const double point[3],
    int quadrature_order,
    double electric_field[3]
);

#endif
