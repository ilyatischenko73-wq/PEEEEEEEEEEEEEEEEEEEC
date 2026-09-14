#ifndef PEEC_RESISTANCE_H
#define PEEC_RESISTANCE_H

#include "dual_mesh.h"
#include "mesh.h"

#include <stddef.h>


/*
 * ============================================================
 * МОДЕЛЬ ПРОВОДНИКА
 * ============================================================
 */
typedef enum
{
    /*
     * Идеальный электрический проводник:
     *
     *     R_s = 0
     *     R_i = 0.
     */
    PEEC_CONDUCTOR_PEC = 0,

    /*
     * Поверхностный резистивный материал.
     *
     * Задается поверхностное сопротивление:
     *
     *     R_s [Ohm].
     *
     * Его также часто называют sheet resistance:
     *
     *     Ohm/square.
     *
     * Численно единица "Ohm/square" имеет размерность Ohm.
     */
    PEEC_CONDUCTOR_SHEET_RESISTANCE

} PeecConductorModel;


/*
 * ============================================================
 * ПАРАМЕТРЫ СОПРОТИВЛЕНИЯ
 * ============================================================
 */
typedef struct
{
    /*
     * Выбранная модель проводника.
     */
    PeecConductorModel model;

    /*
     * Поверхностное сопротивление:
     *
     *     R_s [Ohm].
     *
     * Используется только для:
     *
     *     PEEC_CONDUCTOR_SHEET_RESISTANCE.
     *
     * Для PEC оно игнорируется.
     */
    double sheet_resistance;

    /*
     * Количество OpenMP-потоков.
     */
    int parallel_threads;

} ResistanceOptions;


/*
 * ============================================================
 * МАТРИЦА СОПРОТИВЛЕНИЙ
 * ============================================================
 *
 * В текущей локальной поверхностной модели матрица R
 * диагональна:
 *
 *     R_ij = R_i delta_ij.
 *
 * Поэтому хранить Ne x Ne элементов нет необходимости.
 *
 * Сохраняем только диагональ:
 *
 *     diagonal[i] = R_i [Ohm].
 */
typedef struct
{
    size_t n;

    double *diagonal;

} ResistanceMatrix;


/*
 * Параметры по умолчанию.
 *
 * По умолчанию используем идеальный проводник.
 */
void resistance_options_default(
    ResistanceOptions *options,
    int parallel_threads
);


/*
 * Инициализация.
 */
void resistance_init(
    ResistanceMatrix *R
);


/*
 * Освобождение памяти.
 */
void resistance_free(
    ResistanceMatrix *R
);


/*
 * ============================================================
 * РАСЧЕТ МАТРИЦЫ R
 * ============================================================
 *
 * Для поверхностного сопротивления:
 *
 *     E_t = R_s J_s.
 *
 * На ребре:
 *
 *     J_i = I_i / w_i.
 *
 * Падение напряжения:
 *
 *     U_i^R = R_s * l_i / w_i * I_i.
 *
 * Следовательно:
 *
 *     R_i = R_s * l_i / w_i.
 *
 * Так как:
 *
 *     w_i = S_i^e / l_i,
 *
 * можно записать:
 *
 *     R_i = R_s * l_i^2 / S_i^e.
 *
 *
 * Для PEC:
 *
 *     R_i = 0.
 */
int resistance_compute_matrix(
    const Mesh *mesh,
    const DualMesh *dual,
    const ResistanceOptions *options,
    ResistanceMatrix *R
);


/*
 * Возвращает сопротивление отдельного ребра:
 *
 *     R_i [Ohm].
 */
double resistance_get(
    const ResistanceMatrix *R,
    size_t edge
);


/*
 * Диагностика.
 */
void resistance_print_info(
    const ResistanceMatrix *R,
    const ResistanceOptions *options
);


#endif
