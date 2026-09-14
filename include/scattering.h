#ifndef PEEC_SCATTERING_H
#define PEEC_SCATTERING_H

#include "complex_matrix.h"
#include "incident_field.h"
#include "mesh.h"

#include <stddef.h>


/*
 * ============================================================
 * СПОСОБ ИНТЕГРИРОВАНИЯ ПАДАЮЩЕГО ПОЛЯ ПО РЕБРУ
 * ============================================================
 */
typedef enum
{
    /*
     * Приближение по центру ребра:
     *
     *     U_e ~= E(r_c) . l_e.
     *
     * Работает как для гармонической, так и для
     * временной задачи.
     */
    SCATTERING_EDGE_MIDPOINT = 0,

    /*
     * Точный аналитический интеграл гармонической
     * плоской волны вдоль прямого ребра:
     *
     *     U_e =
     *
     *     [E_hat(r_c) . l_e]
     *
     *     * sinc(
     *           k0 * k_hat.l_e / 2
     *       ).
     *
     *
     * Используется только для гармонической задачи.
     */
    SCATTERING_EDGE_ANALYTIC_PLANE_WAVE = 1

} ScatteringEdgeIntegration;


/*
 * ============================================================
 * ПАРАМЕТРЫ ВОЗБУЖДЕНИЯ
 * ============================================================
 */
typedef struct
{
    ScatteringEdgeIntegration integration;

    int parallel_threads;

} ScatteringOptions;


/*
 * ============================================================
 * ПАРАМЕТРЫ ПО УМОЛЧАНИЮ
 * ============================================================
 *
 * По умолчанию оставляем MIDPOINT, чтобы временная задача
 * продолжала работать без изменения поведения.
 *
 * Для гармонического рассеяния в main.c явно зададим:
 *
 *     SCATTERING_EDGE_ANALYTIC_PLANE_WAVE.
 */
void scattering_options_default(
    ScatteringOptions *options,
    int parallel_threads
);


/*
 * ============================================================
 * ГАРМОНИЧЕСКОЕ ВОЗБУЖДЕНИЕ
 * ============================================================
 *
 * Строит:
 *
 *     U_inc[e] [V].
 *
 *
 * В PEEC-системе:
 *
 *     (R+j*w*L)I + AV = -U_inc.
 *
 *
 * Знак минус здесь НЕ добавляется.
 */
int scattering_build_harmonic_excitation(
    const Mesh *mesh,
    const IncidentField *field,
    const ScatteringOptions *options,
    Complex *edge_voltage
);


/*
 * ============================================================
 * ВРЕМЕННОЕ ВОЗБУЖДЕНИЕ
 * ============================================================
 *
 * Сейчас временная задача использует:
 *
 *     U_e(t) ~= E(r_c,t) . l_e.
 *
 *
 * Аналитическая формула sinc относится к монохроматическому
 * фазору и здесь не применяется.
 */
int scattering_build_time_excitation(
    const Mesh *mesh,
    const IncidentField *field,
    const ScatteringOptions *options,
    double time,
    double *edge_voltage
);


/*
 * ============================================================
 * НУЛЕВОЙ УЗЛОВОЙ ИСТОЧНИК
 * ============================================================
 */
void scattering_zero_node_current_real(
    size_t n_nodes,
    double *node_current
);


void scattering_zero_node_current_complex(
    size_t n_nodes,
    Complex *node_current
);


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void scattering_print_harmonic_excitation_info(
    const Complex *edge_voltage,
    size_t n_edges
);


void scattering_print_time_excitation_info(
    const double *edge_voltage,
    size_t n_edges,
    double time
);


#endif