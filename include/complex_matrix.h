#ifndef PEEC_COMPLEX_MATRIX_H
#define PEEC_COMPLEX_MATRIX_H

#include <stddef.h>


/*
 * ============================================================
 * КОМПЛЕКСНОЕ ЧИСЛО
 * ============================================================
 *
 *     z = re + j*im.
 *
 * Собственный тип используем вместо <complex.h>, чтобы:
 *
 * - полностью контролировать представление данных;
 * - удобно работать с OpenMP;
 * - использовать тот же тип в complex LU;
 * - явно видеть действительную и мнимую части.
 */
typedef struct
{
    double re;
    double im;

} Complex;


/*
 * ============================================================
 * ПЛОТНАЯ КОМПЛЕКСНАЯ МАТРИЦА
 * ============================================================
 *
 * Матрица хранится в row-major формате:
 *
 *     data[row * cols + col].
 *
 * Каждый элемент:
 *
 *     Complex {re, im}.
 */
typedef struct
{
    size_t rows;
    size_t cols;

    Complex *data;

} ComplexMatrix;


/*
 * ============================================================
 * БАЗОВАЯ КОМПЛЕКСНАЯ АРИФМЕТИКА
 * ============================================================
 */

Complex complex_make(
    double re,
    double im
);


Complex complex_add(
    Complex a,
    Complex b
);


Complex complex_sub(
    Complex a,
    Complex b
);


Complex complex_mul(
    Complex a,
    Complex b
);


Complex complex_div(
    Complex a,
    Complex b
);


Complex complex_scale(
    Complex a,
    double scalar
);


Complex complex_conj(
    Complex a
);


double complex_abs(
    Complex a
);


double complex_abs_squared(
    Complex a
);


/*
 * ============================================================
 * МАТРИЦА
 * ============================================================
 */

void complex_matrix_init(
    ComplexMatrix *M
);


int complex_matrix_alloc(
    ComplexMatrix *M,
    size_t rows,
    size_t cols
);


void complex_matrix_free(
    ComplexMatrix *M
);


void complex_matrix_zero(
    ComplexMatrix *M
);


void complex_matrix_fill(
    ComplexMatrix *M,
    Complex value
);


int complex_matrix_copy(
    ComplexMatrix *dst,
    const ComplexMatrix *src
);


Complex complex_matrix_get(
    const ComplexMatrix *M,
    size_t row,
    size_t col
);


void complex_matrix_set(
    ComplexMatrix *M,
    size_t row,
    size_t col,
    Complex value
);


/*
 * ============================================================
 * ТРАНСПОНИРОВАНИЕ
 * ============================================================
 *
 * Обычное транспонирование:
 *
 *     AT = A^T.
 */
int complex_matrix_transpose(
    const ComplexMatrix *A,
    ComplexMatrix *AT,
    int parallel_threads
);


/*
 * ============================================================
 * ЭРМИТОВО СОПРЯЖЕНИЕ
 * ============================================================
 *
 *     AH = A^H = conj(A^T).
 */
int complex_matrix_conjugate_transpose(
    const ComplexMatrix *A,
    ComplexMatrix *AH,
    int parallel_threads
);


/*
 * ============================================================
 * МАТРИЦА x ВЕКТОР
 * ============================================================
 *
 *     y = A*x.
 */
int complex_matrix_vector_multiply(
    const ComplexMatrix *A,
    const Complex *x,
    Complex *y,
    int parallel_threads
);


/*
 * ============================================================
 * МАТРИЦА x МАТРИЦА
 * ============================================================
 *
 *     C = A*B.
 */
int complex_matrix_multiply(
    const ComplexMatrix *A,
    const ComplexMatrix *B,
    ComplexMatrix *C,
    int parallel_threads
);


/*
 * ============================================================
 * ЛИНЕЙНАЯ КОМБИНАЦИЯ
 * ============================================================
 *
 *     C = alpha*A + beta*B.
 */
int complex_matrix_linear_combination(
    const ComplexMatrix *A,
    Complex alpha,
    const ComplexMatrix *B,
    Complex beta,
    ComplexMatrix *C,
    int parallel_threads
);


/*
 * ============================================================
 * НОРМЫ И ДИАГНОСТИКА
 * ============================================================
 */

double complex_matrix_max_abs(
    const ComplexMatrix *A,
    int parallel_threads
);


double complex_matrix_norm_inf(
    const ComplexMatrix *A,
    int parallel_threads
);


/*
 * Для квадратной матрицы:
 *
 *     max |A_ij - A_ji|.
 *
 * Это проверка обычной симметрии:
 *
 *     A = A^T.
 *
 * Не Hermitian symmetry.
 */
double complex_matrix_symmetry_error(
    const ComplexMatrix *A,
    int parallel_threads
);


/*
 * Для квадратной матрицы:
 *
 *     max |A_ij - conj(A_ji)|.
 *
 * Проверка:
 *
 *     A = A^H.
 */
double complex_matrix_hermitian_error(
    const ComplexMatrix *A,
    int parallel_threads
);


void complex_matrix_print_info(
    const ComplexMatrix *M,
    const char *name
);


#endif
