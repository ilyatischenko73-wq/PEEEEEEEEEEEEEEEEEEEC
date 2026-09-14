#include "matrix.h"

#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * ============================================================
 * ПРОВЕРКА РАЗМЕРА ВЫДЕЛЯЕМОЙ ПАМЯТИ
 * ============================================================
 */
static int matrix_check_allocation_size(
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

    if (count > SIZE_MAX / sizeof(double))
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void matrix_init(Matrix *M)
{
    if (M == NULL)
    {
        return;
    }

    memset(M, 0, sizeof(*M));
}


/*
 * ============================================================
 * ВЫДЕЛЕНИЕ ПАМЯТИ
 * ============================================================
 */
int matrix_alloc(
    Matrix *M,
    size_t rows,
    size_t cols
)
{
    if (M == NULL)
    {
        return -1;
    }

    if (matrix_check_allocation_size(rows, cols) != 0)
    {
        fprintf(
            stderr,
            "ERROR: invalid matrix allocation size %zu x %zu.\n",
            rows,
            cols
        );

        return -1;
    }

    matrix_free(M);

    M->rows = rows;
    M->cols = cols;

    M->data = calloc(
        rows * cols,
        sizeof(double)
    );

    if (M->data == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate matrix %zu x %zu.\n",
            rows,
            cols
        );

        matrix_init(M);

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ ПАМЯТИ
 * ============================================================
 */
void matrix_free(Matrix *M)
{
    if (M == NULL)
    {
        return;
    }

    free(M->data);

    matrix_init(M);
}


/*
 * ============================================================
 * ЗАПОЛНЕНИЕ НУЛЯМИ
 * ============================================================
 */
void matrix_zero(Matrix *M)
{
    if (M == NULL || M->data == NULL)
    {
        return;
    }

    memset(
        M->data,
        0,
        M->rows * M->cols * sizeof(double)
    );
}


/*
 * ============================================================
 * ЗАПОЛНЕНИЕ ЗНАЧЕНИЕМ
 * ============================================================
 */
void matrix_fill(
    Matrix *M,
    double value
)
{
    if (M == NULL || M->data == NULL)
    {
        return;
    }

    size_t total = M->rows * M->cols;

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
int matrix_copy(
    Matrix *dst,
    const Matrix *src
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

    /*
     * Копирование объекта в самого себя ничего не меняет.
     */
    if (dst == src)
    {
        return 0;
    }

    if (matrix_alloc(
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
        src->rows * src->cols * sizeof(double)
    );

    return 0;
}


/*
 * ============================================================
 * ДОСТУП К ЭЛЕМЕНТАМ
 * ============================================================
 */
double matrix_get(
    const Matrix *M,
    size_t row,
    size_t col
)
{
    return M->data[row * M->cols + col];
}


void matrix_set(
    Matrix *M,
    size_t row,
    size_t col,
    double value
)
{
    M->data[row * M->cols + col] = value;
}


/*
 * ============================================================
 * ТРАНСПОНИРОВАНИЕ
 * ============================================================
 */
int matrix_transpose(
    const Matrix *A,
    Matrix *AT,
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

    if (matrix_alloc(
        AT,
        A->cols,
        A->rows
    ) != 0)
    {
        return -1;
    }

    /*
     * Каждая строка A независима.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)A->rows;
    ++i_signed)
         {
             size_t i = (size_t)i_signed;

             for (size_t j = 0; j < A->cols; ++j)
             {
                 AT->data[j * AT->cols + i] =
                 A->data[i * A->cols + j];
             }
         }

         return 0;
}


/*
 * ============================================================
 * МАТРИЦА x ВЕКТОР
 * ============================================================
 *
 *     y_i = sum_j A_ij x_j.
 */
int matrix_vector_multiply(
    const Matrix *A,
    const double *x,
    double *y,
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
             size_t i = (size_t)i_signed;

             double sum = 0.0;

             #pragma omp simd reduction(+:sum)
             for (size_t j = 0; j < A->cols; ++j)
             {
                 sum +=
                 A->data[i * A->cols + j]
                 * x[j];
             }

             y[i] = sum;
         }

         return 0;
}


/*
 * ============================================================
 * МАТРИЦА x МАТРИЦА
 * ============================================================
 *
 *     C = A B.
 *
 *
 * Наивный порядок циклов:
 *
 *     i-j-k
 *
 * обычно плохо работает с row-major B, потому что по B
 * происходит скачок между строками.
 *
 * Поэтому используем:
 *
 *     i-k-j.
 *
 * Тогда:
 *
 *     A[i,k]
 *
 * загружается один раз, а B[k,j] читается последовательно
 * по памяти.
 *
 * Для row-major это намного лучше.
 */
int matrix_multiply(
    const Matrix *A,
    const Matrix *B,
    Matrix *C,
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
            "ERROR: incompatible matrix sizes in multiplication: "
            "%zu x %zu and %zu x %zu.\n",
            A->rows,
            A->cols,
            B->rows,
            B->cols
        );

        return -1;
    }

    if (matrix_alloc(
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


    /*
     * Каждая строка C принадлежит только одному потоку.
     *
     * Поэтому нет гонки данных.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)m;
    ++i_signed)
         {
             size_t i = (size_t)i_signed;

             double *c_row =
             &C->data[i * n];

             for (size_t k = 0; k < k_size; ++k)
             {
                 double a_ik =
                 A->data[i * k_size + k];

                 const double *b_row =
                 &B->data[k * n];

                 /*
                  * Если A_ik равен нулю, эта строковая операция
                  * ничего не добавляет.
                  *
                  * Для обычной плотной матрицы это не обязательно
                  * полезно, но для некоторых PEEC-операторов
                  * позволяет пропустить работу.
                  */
                 if (a_ik == 0.0)
                 {
                     continue;
                 }

                 #pragma omp simd
                 for (size_t j = 0; j < n; ++j)
                 {
                     c_row[j] +=
                     a_ik * b_row[j];
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
int matrix_linear_combination(
    const Matrix *A,
    double alpha,
    const Matrix *B,
    double beta,
    Matrix *C,
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

    if (matrix_alloc(
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
             size_t k = (size_t)k_signed;

             C->data[k] =
             alpha * A->data[k]
             + beta * B->data[k];
         }

         return 0;
}


/*
 * ============================================================
 * MATRIX INFINITY NORM
 * ============================================================
 *
 *     ||A||_inf = max_i sum_j |A_ij|.
 */
double matrix_norm_inf(
    const Matrix *A,
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
             size_t i = (size_t)i_signed;

             double row_sum = 0.0;

             #pragma omp simd reduction(+:row_sum)
             for (size_t j = 0; j < A->cols; ++j)
             {
                 row_sum +=
                 fabs(A->data[i * A->cols + j]);
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
 * MAXIMUM ABSOLUTE ELEMENT
 * ============================================================
 */
double matrix_max_abs(
    const Matrix *A,
    int parallel_threads
)
{
    if (A == NULL ||
        A->data == NULL ||
        parallel_threads < 1)
    {
        return NAN;
    }

    double max_value = 0.0;

    size_t total =
    A->rows * A->cols;

    #pragma omp parallel for reduction(max:max_value) schedule(static) num_threads(parallel_threads)
    for (long long k_signed = 0;
         k_signed < (long long)total;
    ++k_signed)
         {
             double value =
             fabs(A->data[(size_t)k_signed]);

             if (value > max_value)
             {
                 max_value = value;
             }
         }

         return max_value;
}


/*
 * ============================================================
 * ПРОВЕРКА СИММЕТРИИ
 * ============================================================
 */
double matrix_symmetry_error(
    const Matrix *A,
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

    size_t n = A->rows;

    double max_error = 0.0;

    #pragma omp parallel for reduction(max:max_error) schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)n;
    ++i_signed)
         {
             size_t i = (size_t)i_signed;

             for (size_t j = i + 1; j < n; ++j)
             {
                 double error = fabs(
                     A->data[i * n + j]
                     - A->data[j * n + i]
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
 * ДИАГНОСТИКА
 * ============================================================
 */
void matrix_print_info(
    const Matrix *M,
    const char *name
)
{
    if (M == NULL)
    {
        return;
    }

    const char *matrix_name =
    (name != NULL) ? name : "MATRIX";

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
        * sizeof(double))
        / (1024.0 * 1024.0);

        printf(
            "Memory : %.3f MiB\n",
            memory_mib
        );
    }

    printf("========================================\n");
}
