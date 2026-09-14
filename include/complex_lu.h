#ifndef PEEC_COMPLEX_LU_H
#define PEEC_COMPLEX_LU_H

#include "complex_matrix.h"

#include <stddef.h>


/*
 * ============================================================
 * ПАРАМЕТРЫ КОМПЛЕКСНОГО LU
 * ============================================================
 */
typedef struct
{
    int parallel_threads;

    /*
     * Относительный порог для определения слишком малого pivot:
     *
     *     |pivot| <= pivot_tolerance * matrix_scale.
     * Масштаб относится к уравновешенной матрице, а не к SI-блокам.
     */
    double pivot_tolerance;

    /*
     * Минимальный размер хвостовой подматрицы,
     * начиная с которого включается OpenMP.
     */
    size_t parallel_threshold;

} ComplexLUOptions;


/*
 * ============================================================
 * КОМПЛЕКСНОЕ LU-РАЗЛОЖЕНИЕ
 * ============================================================
 *
 * Строится:
 *
 *     P A = L U.
 *
 * L и U хранятся совместно:
 *
 *             U U U U
 *             L U U U
 *     data =  L L U U
 *             L L L U
 *
 * Диагональ L равна единице и отдельно не хранится.
 */
typedef struct
{
    size_t n;

    Complex *data;

    size_t *pivots;

    size_t swap_count;

    double matrix_scale;

    /* Aeq = diag(1/row_scale) A diag(1/column_scale). */
    double *row_scale;
    double *column_scale;

    int factorized;

} ComplexLUFactorization;


/*
 * Параметры по умолчанию.
 */
void complex_lu_options_default(
    ComplexLUOptions *options,
    int parallel_threads
);


/*
 * Инициализация.
 */
void complex_lu_init(
    ComplexLUFactorization *lu
);


/*
 * Освобождение памяти.
 */
void complex_lu_free(
    ComplexLUFactorization *lu
);


/*
 * ============================================================
 * ФАКТОРИЗАЦИЯ
 * ============================================================
 *
 * input:
 *
 *     квадратная комплексная матрица n x n.
 *
 * Исходная матрица не изменяется.
 *
 * Возвращает:
 *
 *      0 - успешно;
 *     -1 - ошибка аргументов/памяти;
 *     -2 - матрица численно вырождена.
 */
int complex_lu_factorize(
    const Complex *input,
    size_t n,
    const ComplexLUOptions *options,
    ComplexLUFactorization *lu
);


/*
 * ============================================================
 * РЕШЕНИЕ:
 *
 *     A x = b
 * ============================================================
 */
int complex_lu_solve(
    const ComplexLUFactorization *lu,
    const Complex *b,
    Complex *x
);


/*
 * ============================================================
 * НЕСКОЛЬКО ПРАВЫХ ЧАСТЕЙ
 * ============================================================
 *
 * RHS хранятся последовательно:
 *
 *     B[rhs*n + i].
 */
int complex_lu_solve_multiple(
    const ComplexLUFactorization *lu,
    const Complex *B,
    Complex *X,
    size_t rhs_count,
    int parallel_threads
);


/*
 * Диагностика.
 */
void complex_lu_print_info(
    const ComplexLUFactorization *lu
);


#endif
