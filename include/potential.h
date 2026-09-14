#ifndef PEEC_POTENTIAL_H
#define PEEC_POTENTIAL_H

#include "dual_mesh.h"
#include "mesh.h"

#include <stddef.h>


/*
 * ============================================================
 * ПАРАМЕТРЫ РАСЧЕТА МАТРИЦЫ КОЭФФИЦИЕНТОВ ПОТЕНЦИАЛА
 * ============================================================
 */
typedef struct
{
    /*
     * Порядок квадратуры для хорошо разделенных областей.
     */
    int far_order;

    /*
     * Порядок квадратуры для близких областей.
     */
    int near_order;

    /*
     * Порядок квадратуры для соприкасающихся областей.
     */
    int touching_order;

    /*
     * Порядок SELF-квадратуры.
     */
    int self_order;

    /*
     * Критерий NEAR:
     *
     *     d_min < near_factor * max(h_a, h_b).
     */
    double near_factor;

    /*
     * Количество потоков OpenMP.
     */
    int parallel_threads;

} PotentialOptions;


/*
 * ============================================================
 * МАТРИЦА КОЭФФИЦИЕНТОВ ПОТЕНЦИАЛА
 * ============================================================
 *
 * Размер:
 *
 *     Nv x Nv.
 *
 * Хранение:
 *
 *     row-major.
 *
 * То есть:
 *
 *     P(i,j) = data[i*n+j].
 *
 * Размерность:
 *
 *     В/Кл = 1/Ф.
 */
typedef struct
{
    size_t n;
    double *data;

} PotentialMatrix;


/*
 * Устанавливает параметры по умолчанию.
 */
void potential_options_default(
    PotentialOptions *options,
    int parallel_threads
);


/*
 * Инициализация пустой матрицы.
 */
void potential_init(PotentialMatrix *P);


/*
 * Освобождение памяти.
 */
void potential_free(PotentialMatrix *P);


/*
 * ============================================================
 * РАСЧЕТ ОДНОГО ЭЛЕМЕНТА P_ij
 * ============================================================
 *
 * Используется:
 *
 *              1
 *     P_ij = --------- *
 *            4*pi*eps0
 *
 *              1
 *          -----------
 *          S_i^v S_j^v
 *
 *              /        /
 *             |        |       1
 *             |        |    ------- dS dS'.
 *             |        |    |r-r'|
 *            /Pi_i    /Pi_j
 *
 * Возвращает:
 *
 *      0 - успешно;
 *     -1 - ошибка.
 */
int potential_compute_pair(
    const Mesh *mesh,
    const DualMesh *dual,
    size_t node_i,
    size_t node_j,
    const PotentialOptions *options,
    double *value
);


/*
 * ============================================================
 * РАСЧЕТ ПОЛНОЙ МАТРИЦЫ P
 * ============================================================
 *
 * Матрица симметрична:
 *
 *     P_ij = P_ji.
 *
 * Поэтому вычисляется только верхний треугольник.
 *
 * Расчет распараллелен через OpenMP.
 */
int potential_compute_matrix(
    const Mesh *mesh,
    const DualMesh *dual,
    const PotentialOptions *options,
    PotentialMatrix *P
);


/*
 * Печатает диагностическую информацию.
 */
void potential_print_info(const PotentialMatrix *P);



/*
 * ============================================================
 * СОБСТВЕННЫЕ УЗЛОВЫЕ ЕМКОСТИ F
 * ============================================================
 *
 * Для новой transient-записи:
 *
 *     F_j = 1 / P_jj.
 *
 *
 * F имеет размерность:
 *
 *     [Ф].
 *
 *
 * Вектор F должен быть заранее выделен вызывающей стороной:
 *
 *     double F[P->n].
 */
int potential_build_self_capacitance(
    const PotentialMatrix *P,
    double *F
);


/*
 * ============================================================
 * ПРЕОБРАЗОВАНИЕ V_c -> q
 * ============================================================
 *
 * Собственное узловое напряжение:
 *
 *     V_cj = q_j / F_j.
 *
 * Поэтому:
 *
 *     q_j = F_j V_cj.
 */
int potential_self_voltage_to_charge(
    const PotentialMatrix *P,
    const double *V_c,
    double *charge
);


/*
 * ============================================================
 * ПРЕОБРАЗОВАНИЕ V_c -> phi
 * ============================================================
 *
 * Полный узловой потенциал:
 *
 *     phi = S V_c,
 *
 * где:
 *
 *     S_ja = P_ja / P_aa.
 *
 *
 * Эквивалентно:
 *
 *     q_a   = V_ca / P_aa,
 *
 *     phi_j = sum_a P_ja q_a.
 *
 *
 * Полная матрица S явно не строится.
 */
int potential_self_voltage_to_full_potential(
    const PotentialMatrix *P,
    const double *V_c,
    double *phi
);


#endif
