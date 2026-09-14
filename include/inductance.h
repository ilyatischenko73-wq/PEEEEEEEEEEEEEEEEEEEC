#ifndef PEEC_INDUCTANCE_H
#define PEEC_INDUCTANCE_H

#include "dual_mesh.h"
#include "mesh.h"

#include <stddef.h>


/*
 * ============================================================
 * ПАРАМЕТРЫ РАСЧЕТА МАТРИЦЫ ЧАСТИЧНЫХ ИНДУКТИВНОСТЕЙ
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

    /*
     * Допуск для проверки ортогональности направлений:
     *
     *     |e_a . e_b| < direction_tolerance
     *
     * => L_ab = 0.
     */
    double direction_tolerance;

} InductanceOptions;


/*
 * ============================================================
 * МАТРИЦА ЧАСТИЧНЫХ ИНДУКТИВНОСТЕЙ
 * ============================================================
 *
 * Размер:
 *
 *     Ne x Ne.
 *
 * Хранение:
 *
 *     row-major.
 *
 * То есть:
 *
 *     L(a,b) = data[a*n+b].
 *
 * Размерность:
 *
 *     Гн.
 */
typedef struct
{
    size_t n;
    double *data;

} InductanceMatrix;


/*
 * ============================================================
 * ПАРАМЕТРЫ ПО УМОЛЧАНИЮ
 * ============================================================
 */
void inductance_options_default(
    InductanceOptions *options,
    int parallel_threads
);


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ / ОСВОБОЖДЕНИЕ
 * ============================================================
 */
void inductance_init(
    InductanceMatrix *L
);

void inductance_free(
    InductanceMatrix *L
);


/*
 * ============================================================
 * РАСЧЕТ ОДНОГО КОЭФФИЦИЕНТА L_ab
 * ============================================================
 *
 * Используется:
 *
 *               mu0      e_a . e_b
 *     L_ab = --------- * ----------- *
 *              4*pi        w_a w_b
 *
 *               /        /
 *              |        |       1
 *              |        |    ------- dS dS'.
 *              |        |    |r-r'|
 *             /Pi_a    /Pi_b
 *
 * где:
 *
 *     w_a = S_a^e / l_a,
 *     w_b = S_b^e / l_b.
 *
 * Возвращает:
 *
 *      0 - успешно;
 *     -1 - ошибка.
 */
int inductance_compute_pair(
    const Mesh *mesh,
    const DualMesh *dual,
    size_t edge_a,
    size_t edge_b,
    const InductanceOptions *options,
    double *value
);


/*
 * ============================================================
 * РАСЧЕТ ПОЛНОЙ МАТРИЦЫ L
 * ============================================================
 *
 * Матрица симметрична:
 *
 *     L_ab = L_ba.
 *
 * Поэтому вычисляется только верхний треугольник.
 *
 * Расчет распараллелен через OpenMP.
 */
int inductance_compute_matrix(
    const Mesh *mesh,
    const DualMesh *dual,
    const InductanceOptions *options,
    InductanceMatrix *L
);


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void inductance_print_info(
    const InductanceMatrix *L
);


#endif
