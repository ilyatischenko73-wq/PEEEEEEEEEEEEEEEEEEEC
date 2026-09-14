#include "complex_matrix.h"

#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * ============================================================
 * ВНУТРЕННЯЯ ПРОВЕРКА РАЗМЕРА
 * ============================================================
 */
static int complex_matrix_check_allocation_size(
    size_t rows,
    size_t cols
)
{
    if (rows == 0 || cols == 0)
    {
        return -1;
    }

    if (rows > SIZE_MAX / cols)
    {
        return -1;
    }

    size_t count = rows * cols;

    if (count > SIZE_MAX / sizeof(Complex))
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * КОМПЛЕКСНАЯ АРИФМЕТИКА
 * ============================================================
 */

Complex complex_make(
    double re,
    double im
)
{
    Complex z;

    z.re = re;
    z.im = im;

    return z;
}


Complex complex_add(
    Complex a,
    Complex b
)
{
    Complex z;

    z.re = a.re + b.re;
    z.im = a.im + b.im;

    return z;
}


Complex complex_sub(
    Complex a,
    Complex b
)
{
    Complex z;

    z.re = a.re - b.re;
    z.im = a.im - b.im;

    return z;
}


Complex complex_mul(
    Complex a,
    Complex b
)
{
    /*
     * (a_re + j*a_im)(b_re + j*b_im)
     *
     * =
     *
     * (a_re*b_re - a_im*b_im)
     *
     * + j(a_re*b_im + a_im*b_re).
     */
    Complex z;

    z.re =
    a.re * b.re
    - a.im * b.im;

    z.im =
    a.re * b.im
    + a.im * b.re;

    return z;
}


Complex complex_div(
    Complex a,
    Complex b
)
{
    /*
     *                a * conj(b)
     *     a / b = ----------------
     *                   |b|^2
     */
    double denominator =
    b.re * b.re
    + b.im * b.im;

    Complex z;

    z.re =
    (a.re * b.re + a.im * b.im)
    / denominator;

    z.im =
    (a.im * b.re - a.re * b.im)
    / denominator;

    return z;
}


Complex complex_scale(
    Complex a,
    double scalar
)
{
    Complex z;

    z.re = scalar * a.re;
    z.im = scalar * a.im;

    return z;
}


Complex complex_conj(
    Complex a
)
{
    Complex z;

    z.re = a.re;
    z.im = -a.im;

    return z;
}


double complex_abs_squared(
    Complex a
)
{
    return
    a.re * a.re
    + a.im * a.im;
}


double complex_abs(
    Complex a
)
{
    /*
     * hypot обычно устойчивее, чем:
     *
     *     sqrt(re*re + im*im).
     */
    return hypot(
        a.re,
        a.im
    );
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ МАТРИЦЫ
 * ============================================================
 */
void complex_matrix_init(
    ComplexMatrix *M
)
{
    if (M == NULL)
    {
        return;
    }

    memset(
        M,
        0,
        sizeof(*M)
    );
}


/*
 * ============================================================
 * ВЫДЕЛЕНИЕ ПАМЯТИ
 * ============================================================
 */
int complex_matrix_alloc(
    ComplexMatrix *M,
    size_t rows,
    size_t cols
)
{
    if (M == NULL)
    {
        return -1;
    }

    if (complex_matrix_check_allocation_size(
        rows,
        cols
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: invalid complex matrix size %zu x %zu.\n",
            rows,
            cols
        );

        return -1;
    }

    complex_matrix_free(M);

    M->rows = rows;
    M->cols = cols;

    M->data = calloc(
        rows * cols,
        sizeof(Complex)
    );

    if (M->data == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate complex matrix %zu x %zu.\n",
            rows,
            cols
        );

        complex_matrix_init(M);

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ ПАМЯТИ
 * ============================================================
 */
void complex_matrix_free(
    ComplexMatrix *M
)
{
    if (M == NULL)
    {
        return;
    }

    free(M->data);

    complex_matrix_init(M);
}


/*
 * ============================================================
 * ЗАПОЛНЕНИЕ НУЛЯМИ
 * ============================================================
 */
void complex_matrix_zero(
    ComplexMatrix *M
)
{
    if (M == NULL ||
        M->data == NULL)
    {
        return;
    }

    memset(
        M->data,
        0,
        M->rows * M->cols * sizeof(Complex)
    );
}


/*
 * ============================================================
 * ЗАПОЛНЕНИЕ ЗНАЧЕНИЕМ
 * ============================================================
 */
void complex_matrix_fill(
    ComplexMatrix *M,
    Complex value
)
{
    if (M == NULL ||
        M->data == NULL)
    {
        return;
    }

    size_t total =
    M->rows * M->cols;

    #pragma omp simd
    for (size_t i = 0; i < total; ++i)
    {
        M->data[i] = value;
    }
}


/*
 * ============================================================
 * КОПИРОВАНИЕ
 * ============================================================
 */
int complex_matrix_copy(
    ComplexMatrix *dst,
    const ComplexMatrix *src
)
{
    if (dst == NULL ||
        src == NULL ||
        src->data == NULL ||
        src->rows == 0 ||
        src->cols == 0)
    {
        return -1;
    }

    if (dst == src)
    {
        return 0;
    }

    if (complex_matrix_alloc(
        dst,
        src->rows,
        src->cols
    ) != 0)
    {
        return -1;
    }

    memcpy(
        dst->data,
        src->data,
        src->rows * src->cols * sizeof(Complex)
    );

    return 0;
}


/*
 * ============================================================
 * ДОСТУП
 * ============================================================
 */
Complex complex_matrix_get(
    const ComplexMatrix *M,
    size_t row,
    size_t col
)
{
    return M->data[
        row * M->cols + col
    ];
}


void complex_matrix_set(
    ComplexMatrix *M,
    size_t row,
    size_t col,
    Complex value
)
{
    M->data[
        row * M->cols + col
    ] = value;
}


/*
 * ============================================================
 * ТРАНСПОНИРОВАНИЕ
 * ============================================================
 */
int complex_matrix_transpose(
    const ComplexMatrix *A,
    ComplexMatrix *AT,
    int parallel_threads
)
{
    if (A == NULL ||
        AT == NULL ||
        A->data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (complex_matrix_alloc(
        AT,
        A->cols,
        A->rows
    ) != 0)
    {
        return -1;
    }

    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)A->rows;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             for (size_t j = 0;
                  j < A->cols;
             ++j)
                  {
                      AT->data[
                          j * AT->cols + i
                      ] =
                      A->data[
                          i * A->cols + j
                      ];
                  }
         }

         return 0;
}


/*
 * ============================================================
 * ЭРМИТОВО СОПРЯЖЕНИЕ
 * ============================================================
 */
int complex_matrix_conjugate_transpose(
    const ComplexMatrix *A,
    ComplexMatrix *AH,
    int parallel_threads
)
{
    if (A == NULL ||
        AH == NULL ||
        A->data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (complex_matrix_alloc(
        AH,
        A->cols,
        A->rows
    ) != 0)
    {
        return -1;
    }

    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)A->rows;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             for (size_t j = 0;
                  j < A->cols;
             ++j)
                  {
                      Complex value =
                      A->data[
                          i * A->cols + j
                      ];

                      AH->data[
                          j * AH->cols + i
                      ].re = value.re;

                      AH->data[
                          j * AH->cols + i
                      ].im = -value.im;
                  }
         }

         return 0;
}


/*
 * ============================================================
 * МАТРИЦА x ВЕКТОР
 * ============================================================
 *
 *     y_i = sum_j A_ij*x_j.
 */
int complex_matrix_vector_multiply(
    const ComplexMatrix *A,
    const Complex *x,
    Complex *y,
    int parallel_threads
)
{
    if (A == NULL ||
        A->data == NULL ||
        x == NULL ||
        y == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)A->rows;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             double sum_re = 0.0;
             double sum_im = 0.0;

             /*
              * Если:
              *
              *     A = a + jb
              *     x = c + jd,
              *
              * то:
              *
              *     A*x =
              *
              *     (ac-bd)
              *
              *     + j(ad+bc).
              */
             #pragma omp simd reduction(+:sum_re,sum_im)
             for (size_t j = 0;
                  j < A->cols;
             ++j)
                  {
                      Complex a =
                      A->data[
                          i * A->cols + j
                      ];

                      Complex b =
                      x[j];

                      sum_re +=
                      a.re * b.re
                      - a.im * b.im;

                      sum_im +=
                      a.re * b.im
                      + a.im * b.re;
                  }

                  y[i].re = sum_re;
                  y[i].im = sum_im;
         }

         return 0;
}


/*
 * ============================================================
 * МАТРИЦА x МАТРИЦА
 * ============================================================
 *
 *     C = A*B.
 *
 * Используем порядок:
 *
 *     i-k-j,
 *
 * чтобы строка B[k,:] читалась последовательно.
 */
int complex_matrix_multiply(
    const ComplexMatrix *A,
    const ComplexMatrix *B,
    ComplexMatrix *C,
    int parallel_threads
)
{
    if (A == NULL ||
        B == NULL ||
        C == NULL ||
        A->data == NULL ||
        B->data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (A->cols != B->rows)
    {
        fprintf(
            stderr,
            "ERROR: incompatible complex matrix sizes: "
            "%zu x %zu and %zu x %zu.\n",
            A->rows,
            A->cols,
            B->rows,
            B->cols
        );

        return -1;
    }

    if (complex_matrix_alloc(
        C,
        A->rows,
        B->cols
    ) != 0)
    {
        return -1;
    }

    size_t m = A->rows;
    size_t k_size = A->cols;
    size_t n = B->cols;

    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)m;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             Complex *c_row =
             &C->data[i * n];

             for (size_t k = 0;
                  k < k_size;
             ++k)
                  {
                      Complex a =
                      A->data[
                          i * k_size + k
                      ];

                      /*
                       * Если элемент точно нулевой,
                       * эту операцию можно пропустить.
                       */
                      if (a.re == 0.0 &&
                          a.im == 0.0)
                      {
                          continue;
                      }

                      const Complex *b_row =
                      &B->data[k * n];

                      #pragma omp simd
                      for (size_t j = 0;
                           j < n;
                      ++j)
                           {
                               double re =
                               a.re * b_row[j].re
                               - a.im * b_row[j].im;

                               double im =
                               a.re * b_row[j].im
                               + a.im * b_row[j].re;

                               c_row[j].re += re;
                               c_row[j].im += im;
                           }
                  }
         }

         return 0;
}


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
)
{
    if (A == NULL ||
        B == NULL ||
        C == NULL ||
        A->data == NULL ||
        B->data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (A->rows != B->rows ||
        A->cols != B->cols)
    {
        return -1;
    }

    if (complex_matrix_alloc(
        C,
        A->rows,
        A->cols
    ) != 0)
    {
        return -1;
    }

    size_t total =
    A->rows * A->cols;

    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long k_signed = 0;
         k_signed < (long long)total;
    ++k_signed)
         {
             size_t k =
             (size_t)k_signed;

             Complex a =
             A->data[k];

             Complex b =
             B->data[k];

             Complex alpha_a =
             complex_mul(
                 alpha,
                 a
             );

             Complex beta_b =
             complex_mul(
                 beta,
                 b
             );

             C->data[k].re =
             alpha_a.re + beta_b.re;

             C->data[k].im =
             alpha_a.im + beta_b.im;
         }

         return 0;
}


/*
 * ============================================================
 * MAX |A_ij|
 * ============================================================
 */
double complex_matrix_max_abs(
    const ComplexMatrix *A,
    int parallel_threads
)
{
    if (A == NULL ||
        A->data == NULL ||
        parallel_threads < 1)
    {
        return NAN;
    }

    size_t total =
    A->rows * A->cols;

    double max_value = 0.0;

    #pragma omp parallel for reduction(max:max_value) schedule(static) num_threads(parallel_threads)
    for (long long k_signed = 0;
         k_signed < (long long)total;
    ++k_signed)
         {
             Complex value =
             A->data[(size_t)k_signed];

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
 * INFINITY NORM
 * ============================================================
 *
 *     ||A||_inf =
 *
 *     max_i sum_j |A_ij|.
 */
double complex_matrix_norm_inf(
    const ComplexMatrix *A,
    int parallel_threads
)
{
    if (A == NULL ||
        A->data == NULL ||
        parallel_threads < 1)
    {
        return NAN;
    }

    double max_row_sum = 0.0;

    #pragma omp parallel for reduction(max:max_row_sum) schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)A->rows;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             double row_sum = 0.0;

             #pragma omp simd reduction(+:row_sum)
             for (size_t j = 0;
                  j < A->cols;
             ++j)
                  {
                      Complex value =
                      A->data[
                          i * A->cols + j
                      ];

                      row_sum +=
                      hypot(
                          value.re,
                          value.im
                      );
                  }

                  if (row_sum > max_row_sum)
                  {
                      max_row_sum = row_sum;
                  }
         }

         return max_row_sum;
}


/*
 * ============================================================
 * ОШИБКА СИММЕТРИИ
 * ============================================================
 *
 *     max |A_ij - A_ji|.
 */
double complex_matrix_symmetry_error(
    const ComplexMatrix *A,
    int parallel_threads
)
{
    if (A == NULL ||
        A->data == NULL ||
        A->rows != A->cols ||
        parallel_threads < 1)
    {
        return NAN;
    }

    size_t n =
    A->rows;

    double max_error = 0.0;

    #pragma omp parallel for reduction(max:max_error) schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)n;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             for (size_t j = i + 1;
                  j < n;
             ++j)
                  {
                      Complex a =
                      A->data[
                          i * n + j
                      ];

                      Complex b =
                      A->data[
                          j * n + i
                      ];

                      double error_re =
                      a.re - b.re;

                      double error_im =
                      a.im - b.im;

                      double error =
                      hypot(
                          error_re,
                          error_im
                      );

                      if (error > max_error)
                      {
                          max_error = error;
                      }
                  }
         }

         return max_error;
}


/*
 * ============================================================
 * ОШИБКА ЭРМИТОВОСТИ
 * ============================================================
 *
 *     max |A_ij - conj(A_ji)|.
 */
double complex_matrix_hermitian_error(
    const ComplexMatrix *A,
    int parallel_threads
)
{
    if (A == NULL ||
        A->data == NULL ||
        A->rows != A->cols ||
        parallel_threads < 1)
    {
        return NAN;
    }

    size_t n =
    A->rows;

    double max_error = 0.0;

    /*
     * В отличие от проверки обычной симметрии,
     * здесь нужно также проверить диагональ:
     *
     * Hermitian-матрица должна иметь действительную диагональ.
     */
    #pragma omp parallel for reduction(max:max_error) schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)n;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             /*
              * Для диагонали:
              *
              *     A_ii = conj(A_ii)
              *
              * поэтому Im(A_ii) = 0.
              *
              * Разность:
              *
              *     A_ii - conj(A_ii)
              *
              * имеет модуль:
              *
              *     2*|Im(A_ii)|.
              */
             double diagonal_error =
             2.0 * fabs(
                 A->data[i * n + i].im
             );

             if (diagonal_error > max_error)
             {
                 max_error =
                 diagonal_error;
             }

             for (size_t j = i + 1;
                  j < n;
             ++j)
                  {
                      Complex a =
                      A->data[
                          i * n + j
                      ];

                      Complex b =
                      A->data[
                          j * n + i
                      ];

                      /*
                       * a - conj(b):
                       *
                       * real:
                       *
                       *     a.re - b.re
                       *
                       * imag:
                       *
                       *     a.im + b.im.
                       */
                      double error_re =
                      a.re - b.re;

                      double error_im =
                      a.im + b.im;

                      double error =
                      hypot(
                          error_re,
                          error_im
                      );

                      if (error > max_error)
                      {
                          max_error =
                          error;
                      }
                  }
         }

         return max_error;
}


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void complex_matrix_print_info(
    const ComplexMatrix *M,
    const char *name
)
{
    if (M == NULL)
    {
        return;
    }

    const char *matrix_name =
    (name != NULL)
    ? name
    : "COMPLEX MATRIX";

    printf("\n");
    printf("========================================\n");
    printf("%s\n", matrix_name);
    printf("========================================\n");
    printf(
        "Size   : %zu x %zu\n",
        M->rows,
        M->cols
    );

    if (M->data != NULL)
    {
        double memory_mib =
        ((double)M->rows
        * (double)M->cols
        * sizeof(Complex))
        / (1024.0 * 1024.0);

        printf(
            "Memory : %.3f MiB\n",
            memory_mib
        );
    }

    printf("========================================\n");
}
