#include "complex_lu.h"

#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * ============================================================
 * ПАРАМЕТРЫ ПО УМОЛЧАНИЮ
 * ============================================================
 */
void complex_lu_options_default(
    ComplexLUOptions *options,
    int parallel_threads
)
{
    if (options == NULL)
    {
        return;
    }

    options->parallel_threads = parallel_threads;
    options->pivot_tolerance = 1.0e-14;
    options->parallel_threshold = 64;
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void complex_lu_init(
    ComplexLUFactorization *lu
)
{
    if (lu == NULL)
    {
        return;
    }

    memset(
        lu,
        0,
        sizeof(*lu)
    );
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ ПАМЯТИ
 * ============================================================
 */
void complex_lu_free(
    ComplexLUFactorization *lu
)
{
    if (lu == NULL)
    {
        return;
    }

    free(lu->data);
    free(lu->pivots);
    free(lu->row_scale);
    free(lu->column_scale);

    complex_lu_init(lu);
}


/*
 * ============================================================
 * ПЕРЕСТАНОВКА ДВУХ СТРОК
 * ============================================================
 */
static void complex_swap_rows(
    Complex *matrix,
    size_t n,
    size_t row_a,
    size_t row_b
)
{
    if (row_a == row_b)
    {
        return;
    }

    Complex *a =
    &matrix[row_a * n];

    Complex *b =
    &matrix[row_b * n];

    #pragma omp simd
    for (size_t j = 0; j < n; ++j)
    {
        Complex temp = a[j];
        a[j] = b[j];
        b[j] = temp;
    }
}


/*
 * ============================================================
 * МАКСИМАЛЬНЫЙ МОДУЛЬ ЭЛЕМЕНТА МАТРИЦЫ
 * ============================================================
 *
 *     scale = max |A_ij|.
 */
static double complex_matrix_scale(
    const Complex *matrix,
    size_t n,
    int parallel_threads
)
{
    size_t total =
    n * n;

    double max_value = 0.0;

    #pragma omp parallel for reduction(max:max_value) schedule(static) num_threads(parallel_threads) if(total >= 65536)
    for (long long k_signed = 0;
         k_signed < (long long)total;
    ++k_signed)
         {
             Complex value =
             matrix[(size_t)k_signed];

             double magnitude =
             hypot(
                 value.re,
                 value.im
             );

             if (magnitude > max_value)
             {
                 max_value = magnitude;
             }
         }

         return max_value;
}


/*
 * ============================================================
 * КОМПЛЕКСНОЕ ДЕЛЕНИЕ БЕЗ ВЫЗОВА ВНЕШНЕЙ ФУНКЦИИ
 * ============================================================
 *
 *     a / b.
 *
 * Используем прямо внутри LU, чтобы сократить накладные
 * расходы в основном цикле.
 */
static inline Complex complex_div_inline(
    Complex a,
    Complex b
)
{
    double denominator =
    b.re * b.re
    + b.im * b.im;

    Complex z;

    z.re =
    (a.re * b.re
    + a.im * b.im)
    / denominator;

    z.im =
    (a.im * b.re
    - a.re * b.im)
    / denominator;

    return z;
}


/*
 * ============================================================
 * LU-ФАКТОРИЗАЦИЯ
 * ============================================================
 *
 * На шаге k:
 *
 * 1. ищем pivot:
 *
 *        max |A_ik|
 *
 * 2. переставляем строки;
 *
 * 3. вычисляем:
 *
 *        L_ik = A_ik / A_kk
 *
 * 4. обновляем:
 *
 *        A_ij <- A_ij - L_ik*A_kj.
 */
int complex_lu_factorize(
    const Complex *input,
    size_t n,
    const ComplexLUOptions *options,
    ComplexLUFactorization *lu
)
{
    if (input == NULL ||
        options == NULL ||
        lu == NULL ||
        n == 0 ||
        options->parallel_threads < 1 ||
        !isfinite(options->pivot_tolerance) || options->pivot_tolerance < 0.0)
    {
        return -1;
    }

    if (n > SIZE_MAX / n)
    {
        fprintf(
            stderr,
            "ERROR: complex LU matrix size overflow.\n"
        );

        return -1;
    }

    size_t element_count =
    n * n;

    if (element_count > SIZE_MAX / sizeof(Complex))
    {
        fprintf(
            stderr,
            "ERROR: complex LU allocation size overflow.\n"
        );

        return -1;
    }

    complex_lu_init(lu);

    lu->n = n;

    lu->data = malloc(
        element_count * sizeof(Complex)
    );

    lu->pivots = malloc(
        n * sizeof(size_t)
    );

    lu->row_scale = malloc(n * sizeof(double));
    lu->column_scale = malloc(n * sizeof(double));
    if (lu->data == NULL ||
        lu->pivots == NULL || lu->row_scale == NULL || lu->column_scale == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate complex LU factorization.\n"
        );

        complex_lu_free(lu);

        return -1;
    }

    /* Equilibrate before selecting a relative pivot threshold.  PEEC mixes
     * volts, currents, inductances and potential coefficients in SI units. */
    for (size_t i = 0; i < n; ++i) {
        double scale = 0.0;
        for (size_t j = 0; j < n; ++j) {
            Complex z = input[i*n+j];
            double magnitude = hypot(z.re,z.im);
            if (!isfinite(magnitude)) { complex_lu_free(lu); return -2; }
            scale = fmax(scale,magnitude);
        }
        if (scale == 0.0) { complex_lu_free(lu); return -2; }
        lu->row_scale[i] = scale;
        for (size_t j = 0; j < n; ++j) {
            lu->data[i*n+j].re = input[i*n+j].re / scale;
            lu->data[i*n+j].im = input[i*n+j].im / scale;
        }
    }
    for (size_t j = 0; j < n; ++j) {
        double scale = 0.0;
        for (size_t i = 0; i < n; ++i)
            scale = fmax(scale,hypot(lu->data[i*n+j].re,lu->data[i*n+j].im));
        if (scale == 0.0) { complex_lu_free(lu); return -2; }
        lu->column_scale[j] = scale;
        for (size_t i = 0; i < n; ++i) {
            lu->data[i*n+j].re /= scale;
            lu->data[i*n+j].im /= scale;
        }
    }
    lu->matrix_scale = complex_matrix_scale(lu->data,n,options->parallel_threads);

    if (!isfinite(lu->matrix_scale) ||
        lu->matrix_scale <= 0.0)
    {
        fprintf(
            stderr,
            "ERROR: complex LU input matrix is invalid or zero.\n"
        );

        complex_lu_free(lu);

        return -2;
    }

    double pivot_limit =
    options->pivot_tolerance
    * lu->matrix_scale;

    lu->swap_count = 0;


    /*
     * ========================================================
     * ОСНОВНОЙ ЦИКЛ
     * ========================================================
     */
    for (size_t k = 0; k < n; ++k)
    {
        /*
         * ----------------------------------------------------
         * ПОИСК PIVOT
         * ----------------------------------------------------
         */
        size_t pivot_row = k;

        Complex pivot_candidate =
        lu->data[k * n + k];

        double pivot_abs =
        hypot(
            pivot_candidate.re,
            pivot_candidate.im
        );

        for (size_t i = k + 1; i < n; ++i)
        {
            Complex value =
            lu->data[i * n + k];

            double magnitude =
            hypot(
                value.re,
                value.im
            );

            if (magnitude > pivot_abs)
            {
                pivot_abs = magnitude;
                pivot_row = i;
            }
        }


        /*
         * ----------------------------------------------------
         * ПРОВЕРКА PIVOT
         * ----------------------------------------------------
         */
        if (!isfinite(pivot_abs) ||
            pivot_abs <= pivot_limit)
        {
            fprintf(
                stderr,
                "ERROR: complex LU matrix is singular or nearly "
                "singular at column %zu, |pivot| = %.9e.\n",
                k,
                pivot_abs
            );

            lu->factorized = 0;

            return -2;
        }


        /*
         * Сохраняем историю перестановок.
         */
        lu->pivots[k] =
        pivot_row;


        /*
         * ----------------------------------------------------
         * ПЕРЕСТАНОВКА СТРОК
         * ----------------------------------------------------
         */
        if (pivot_row != k)
        {
            complex_swap_rows(
                lu->data,
                n,
                k,
                pivot_row
            );

            ++lu->swap_count;
        }


        /*
         * Последняя строка не имеет хвостовой подматрицы.
         */
        if (k + 1 >= n)
        {
            continue;
        }


        Complex pivot =
        lu->data[k * n + k];

        double pivot_magnitude =
        hypot(
            pivot.re,
            pivot.im
        );

        if (!isfinite(pivot_magnitude) ||
            pivot_magnitude <= pivot_limit)
        {
            fprintf(
                stderr,
                "ERROR: invalid complex LU pivot at step %zu.\n",
                k
            );

            lu->factorized = 0;

            return -2;
        }


        /*
         * ----------------------------------------------------
         * ELIMINATION + SCHUR UPDATE
         * ----------------------------------------------------
         */
        size_t remaining =
        n - k - 1;

        #pragma omp parallel for schedule(static) num_threads(options->parallel_threads) if(remaining >= options->parallel_threshold)
        for (long long i_signed = (long long)k + 1;
             i_signed < (long long)n;
        ++i_signed)
             {
                 size_t i =
                 (size_t)i_signed;

                 Complex a_ik =
                 lu->data[i * n + k];

                 Complex multiplier =
                 complex_div_inline(
                     a_ik,
                     pivot
                 );

                 /*
                  * Сохраняем:
                  *
                  *     L_ik.
                  */
                 lu->data[i * n + k] =
                 multiplier;


                 /*
                  *     A_ij -= L_ik * U_kj.
                  */
                 #pragma omp simd
                 for (size_t j = k + 1;
                      j < n;
                 ++j)
                      {
                          Complex u =
                          lu->data[k * n + j];

                          /*
                           * multiplier*u.
                           */
                          double product_re =
                          multiplier.re * u.re
                          - multiplier.im * u.im;

                          double product_im =
                          multiplier.re * u.im
                          + multiplier.im * u.re;

                          lu->data[i * n + j].re -=
                          product_re;

                          lu->data[i * n + j].im -=
                          product_im;
                      }
             }
    }


    /*
     * ========================================================
     * ФИНАЛЬНАЯ ПРОВЕРКА ДИАГОНАЛИ U
     * ========================================================
     */
    for (size_t i = 0; i < n; ++i)
    {
        Complex diagonal =
        lu->data[i * n + i];

        double magnitude =
        hypot(
            diagonal.re,
            diagonal.im
        );

        if (!isfinite(magnitude) ||
            magnitude <= pivot_limit)
        {
            fprintf(
                stderr,
                "ERROR: invalid complex U diagonal "
                "at [%zu,%zu], |Uii| = %.9e.\n",
                i,
                i,
                magnitude
            );

            lu->factorized = 0;

            return -2;
        }
    }

    lu->factorized = 1;

    return 0;
}


/*
 * ============================================================
 * ПРИМЕНЕНИЕ ПЕРЕСТАНОВОК
 * ============================================================
 *
 *     b <- P*b.
 */
static void complex_apply_pivots(
    const ComplexLUFactorization *lu,
    Complex *x
)
{
    for (size_t k = 0; k < lu->n; ++k)
    {
        size_t pivot_row =
        lu->pivots[k];

        if (pivot_row != k)
        {
            Complex temp =
            x[k];

            x[k] =
            x[pivot_row];

            x[pivot_row] =
            temp;
        }
    }
}


/*
 * ============================================================
 * FORWARD SUBSTITUTION
 * ============================================================
 *
 * Решаем:
 *
 *     L y = P b.
 *
 * Диагональ L:
 *
 *     L_ii = 1.
 */
static void complex_forward_substitution(
    const ComplexLUFactorization *lu,
    Complex *x
)
{
    size_t n =
    lu->n;

    for (size_t i = 0; i < n; ++i)
    {
        double correction_re = 0.0;
        double correction_im = 0.0;

        #pragma omp simd reduction(+:correction_re,correction_im)
        for (size_t j = 0; j < i; ++j)
        {
            Complex l =
            lu->data[i * n + j];

            Complex value =
            x[j];

            correction_re +=
            l.re * value.re
            - l.im * value.im;

            correction_im +=
            l.re * value.im
            + l.im * value.re;
        }

        x[i].re -= correction_re;
        x[i].im -= correction_im;
    }
}


/*
 * ============================================================
 * BACKWARD SUBSTITUTION
 * ============================================================
 *
 * Решаем:
 *
 *     U x = y.
 */
static int complex_backward_substitution(
    const ComplexLUFactorization *lu,
    Complex *x
)
{
    size_t n =
    lu->n;

    for (size_t ii = n; ii-- > 0;)
    {
        size_t i =
        ii;

        double correction_re = 0.0;
        double correction_im = 0.0;

        #pragma omp simd reduction(+:correction_re,correction_im)
        for (size_t j = i + 1;
             j < n;
        ++j)
             {
                 Complex u =
                 lu->data[i * n + j];

                 Complex value =
                 x[j];

                 correction_re +=
                 u.re * value.re
                 - u.im * value.im;

                 correction_im +=
                 u.re * value.im
                 + u.im * value.re;
             }

             Complex rhs;

             rhs.re =
             x[i].re - correction_re;

             rhs.im =
             x[i].im - correction_im;

             Complex diagonal =
             lu->data[i * n + i];

             double denominator =
             diagonal.re * diagonal.re
             + diagonal.im * diagonal.im;

             if (!isfinite(denominator) ||
                 denominator == 0.0)
             {
                 return -1;
             }

             /*
              * rhs / diagonal.
              */
             x[i].re =
             (rhs.re * diagonal.re
             + rhs.im * diagonal.im)
             / denominator;

             x[i].im =
             (rhs.im * diagonal.re
             - rhs.re * diagonal.im)
             / denominator;
    }

    return 0;
}


/*
 * ============================================================
 * РЕШЕНИЕ ОДНОЙ СИСТЕМЫ
 * ============================================================
 */
int complex_lu_solve(
    const ComplexLUFactorization *lu,
    const Complex *b,
    Complex *x
)
{
    if (lu == NULL ||
        b == NULL ||
        x == NULL ||
        lu->data == NULL ||
        lu->pivots == NULL ||
        !lu->factorized)
    {
        return -1;
    }

    size_t n =
    lu->n;

    for (size_t i = 0; i < n; ++i) {
        x[i].re = b[i].re / lu->row_scale[i];
        x[i].im = b[i].im / lu->row_scale[i];
    }

    /*
     *     PA = LU
     *
     *     Ax = b
     *
     * =>  LUx = Pb.
     */
    complex_apply_pivots(
        lu,
        x
    );

    complex_forward_substitution(
        lu,
        x
    );

    if (complex_backward_substitution(
        lu,
        x
    ) != 0)
    {
        return -1;
    }

    for (size_t i = 0; i < n; ++i) {
        x[i].re /= lu->column_scale[i];
        x[i].im /= lu->column_scale[i];
    }
    return 0;
}


/*
 * ============================================================
 * НЕСКОЛЬКО ПРАВЫХ ЧАСТЕЙ
 * ============================================================
 */
int complex_lu_solve_multiple(
    const ComplexLUFactorization *lu,
    const Complex *B,
    Complex *X,
    size_t rhs_count,
    int parallel_threads
)
{
    if (lu == NULL ||
        B == NULL ||
        X == NULL ||
        rhs_count == 0 ||
        parallel_threads < 1 ||
        !lu->factorized)
    {
        return -1;
    }

    int error_code = 0;

    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long rhs_signed = 0;
         rhs_signed < (long long)rhs_count;
    ++rhs_signed)
         {
             size_t rhs =
             (size_t)rhs_signed;

             const Complex *b =
             &B[rhs * lu->n];

             Complex *x =
             &X[rhs * lu->n];

             int status =
             complex_lu_solve(
                 lu,
                 b,
                 x
             );

             if (status != 0)
             {
                 #pragma omp atomic write
                 error_code = 1;
             }
         }

         if (error_code != 0)
         {
             return -1;
         }

         return 0;
}


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void complex_lu_print_info(
    const ComplexLUFactorization *lu
)
{
    if (lu == NULL)
    {
        return;
    }

    printf("\n");
    printf("========================================\n");
    printf("COMPLEX LU FACTORIZATION\n");
    printf("========================================\n");
    printf(
        "Size         : %zu x %zu\n",
        lu->n,
        lu->n
    );
    printf(
        "Swap count   : %zu\n",
        lu->swap_count
    );
    printf(
        "Matrix scale : %.9e\n",
        lu->matrix_scale
    );
    printf(
        "Factorized   : %s\n",
        lu->factorized
        ? "yes"
        : "no"
    );

    if (lu->data != NULL &&
        lu->n > 0)
    {
        double memory_mib =
        ((double)lu->n
        * (double)lu->n
        * sizeof(Complex))
        / (1024.0 * 1024.0);

        printf(
            "LU memory    : %.3f MiB\n",
            memory_mib
        );
    }

    printf("========================================\n");
}
