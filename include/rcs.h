#ifndef PEEC_RCS_H
#define PEEC_RCS_H

#include "complex_matrix.h"
#include "dual_mesh.h"
#include "harmonic.h"
#include "incident_field.h"
#include "mesh.h"

#include <stddef.h>


/*
 * ============================================================
 * ИНТЕГРИРОВАНИЕ ТОКА ДЛЯ ДАЛЬНЕЙ ЗОНЫ
 * ============================================================
 */
typedef enum
{
    /*
     * Старое midpoint-приближение:
     *
     *     F_e ~= I_e * l_e_vec
     *            * exp(j*k*r_hat.r_e).
     */
    RCS_CURRENT_MIDPOINT = 0,

    /*
     * Интегрирование непосредственно по дуальной
     * поверхности ребра Pi_e^e:
     *
     *     G_e =
     *
     *     integral_{Pi_e^e}
     *     exp(j*k*r_hat.r') dS'.
     *
     * Затем:
     *
     *     F_e =
     *
     *     I_e / w_e
     *     * e_e
     *     * G_e.
     */
    RCS_CURRENT_DUAL_PATCH_QUADRATURE = 1

} RcsIntegration;


/*
 * ============================================================
 * НАСТРОЙКИ RCS
 * ============================================================
 */
typedef struct
{
    RcsIntegration integration;

    /*
     * Порядок Gauss-Legendre для гладкого поверхностного
     * интеграла по каждому треугольнику.
     *
     * Поддерживаются те же значения, что quadrature.c:
     *
     *     2,3,4,5,6,8,10,12,16.
     */
    int surface_order;

    int parallel_threads;

} RcsOptions;


/*
 * ============================================================
 * РЕЗУЛЬТАТ ДЛЯ ОДНОГО НАПРАВЛЕНИЯ
 * ============================================================
 */
typedef struct
{
    double theta_deg;
    double phi_deg;

    /*
     * Направление от объекта к наблюдателю:
     *
     *     r_hat.
     */
    double direction[3];

    /*
     * Radiation integral:
     *
     *     F =
     *
     *     integral J(r')
     *     exp(j*k*r_hat.r')
     *     dS'.
     *
     * [A*m]
     */
    Complex radiation_vector[3];

    /*
     *     T =
     *
     *     r_hat x (r_hat x F).
     */
    Complex transverse_vector[3];

    /*
     * RCS [m^2].
     */
    double sigma;

    /*
     * RCS [dBsm].
     */
    double sigma_dbsm;

} RcsResult;


/*
 * ============================================================
 * НАСТРОЙКИ
 * ============================================================
 */
void rcs_options_default(
    RcsOptions *options,
    int parallel_threads
);


/*
 * ============================================================
 * НАПРАВЛЕНИЕ НАБЛЮДЕНИЯ
 * ============================================================
 */
int rcs_build_observation_direction(
    double theta_deg,
    double phi_deg,
    double direction[3]
);


/*
 * ============================================================
 * RADIATION VECTOR
 * ============================================================
 */
int rcs_compute_radiation_vector(
    const Mesh *mesh,
    const DualMesh *dual,
    const Complex *edge_current,
    double wave_number,
    const double observation_direction[3],
    const RcsOptions *options,
    Complex radiation_vector[3]
);


/*
 * ============================================================
 * RCS В ОДНОМ НАПРАВЛЕНИИ
 * ============================================================
 */
int rcs_compute_direction(
    const Mesh *mesh,
    const DualMesh *dual,
    const IncidentField *incident_field,
    const HarmonicSolution *solution,
    const RcsOptions *options,
    double theta_deg,
    double phi_deg,
    RcsResult *result
);


/*
 * ============================================================
 * ДАЛЬНЕЕ ЭЛЕКТРИЧЕСКОЕ ПОЛЕ
 * ============================================================
 */
int rcs_compute_far_electric_field(
    const RcsResult *result,
    double wave_number,
    double distance,
    Complex E_sca[3]
);


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void rcs_print_result(
    const RcsResult *result
);


#endif