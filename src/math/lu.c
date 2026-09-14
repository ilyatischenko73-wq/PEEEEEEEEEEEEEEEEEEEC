#include "lu.h"

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
void lu_options_default(
    LUOptions *options,
    int parallel_threads
)
{
    if (options == NULL)
    {
        return;
    }

    /*
     * Число OpenMP-потоков приходит из:
     *
     *     --parallel N.
     */
    options->parallel_threads = parallel_threads;

    /*
     * Относительный порог определения численно нулевого pivot.
     *
     * Это не абсолютный порог, а множитель масштаба матрицы.
     */
    options->pivot_tolerance = 1.0e-14;

    /*
     * Для маленькой хвостовой матрицы OpenMP может только
     * замедлить вычисления.
     *
     * 64 - разумная начальная граница, которую потом можно
     * подобрать экспериментально.
     */
    options->parallel_threshold = 64;
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void lu_init(LUFactorization *lu)
{
    if (lu == NULL)
    {
        return;
    }

    memset(lu, 0, sizeof(*lu));
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ ПАМЯТИ
 * ============================================================
 */
void lu_free(LUFactorization *lu)
{
    if (lu == NULL)
    {
        return;
    }

    free(lu->data);
    free(lu->pivots);

    lu_init(lu);
}


/*
 * ============================================================
 * ПЕРЕСТАНОВКА ДВУХ СТРОК
 * ============================================================
 *
 * Матрица хранится в row-major формате.
 */
static void swap_rows(
    double *matrix,
    size_t n,
    size_t row_a,
    size_t row_b
)
{
    if (row_a == row_b)
    {
        return;
    }

    double *a = &matrix[row_a * n];
    double *b = &matrix[row_b * n];

    /*
     * Элементы строк независимы.
     *
     * Здесь OpenMP parallel не нужен: операция слишком дешевая.
     * Но разрешаем компилятору векторизовать цикл.
     */
    #pragma omp simd
    for (size_t j = 0; j < n; ++j)
    {
        double temp = a[j];
        a[j] = b[j];
        b[j] = temp;
    }
}


/*
 * ============================================================
 * МАСШТАБ МАТРИЦЫ
 * ============================================================
 *
 * Вычисляем:
 *
 *     scale = max |A_ij|.
 *
 * Это позволит использовать относительный критерий
 * вырожденности.
 */
static double matrix_max_abs(
    const double *matrix,
    size_t n,
    int parallel_threads
)
{
    double max_value = 0.0;

    size_t total = n * n;

    /*
     * reduction(max:...) поддерживается OpenMP для double.
     */
    #pragma omp parallel for reduction(max:max_value) schedule(static) num_threads(parallel_threads) if(total >= 65536)
    for (long long k = 0; k < (long long)total; ++k)
    {
        double value = fabs(matrix[(size_t)k]);

        if (value > max_value)
        {
            max_value = value;
        }
    }

    return max_value;
}


/*
 * ============================================================
 * LU-ФАКТОРИЗАЦИЯ С ЧАСТИЧНЫМ ВЫБОРОМ ГЛАВНОГО ЭЛЕМЕНТА
 * ============================================================
 *
 * На шаге k:
 *
 * 1. В столбце k ищем:
 *
 *        pivot_row = argmax |A_ik|,
 *
 *        i = k,...,n-1.
 *
 * 2. Меняем строки k и pivot_row.
 *
 * 3. Вычисляем множители:
 *
 *        L_ik = A_ik / A_kk.
 *
 * 4. Обновляем хвостовую матрицу:
 *
 *        A_ij <- A_ij - L_ik A_kj.
 *
 *
 * Именно шаг 4 содержит основную вычислительную работу:
 *
 *        O(n^3).
 *
 * Он распараллелен по строкам i.
 */
int lu_factorize(
    const double *input,
    size_t n,
    const LUOptions *options,
    LUFactorization *lu
)
{
    if (input == NULL ||
        options == NULL ||
        lu == NULL ||
        n == 0)
    {
        return -1;
    }

    if (options->parallel_threads < 1)
    {
        return -1;
    }

    /*
     * Проверяем переполнение при вычислении n*n.
     */
    if (n > SIZE_MAX / n)
    {
        fprintf(stderr, "ERROR: LU matrix size overflow.\n");
        return -1;
    }

    size_t element_count = n * n;

    if (element_count > SIZE_MAX / sizeof(double))
    {
        fprintf(stderr, "ERROR: LU allocation size overflow.\n");
        return -1;
    }

    /*
     * Если структура ранее использовалась, пользователь должен
     * был вызвать lu_free().
     *
     * Здесь начинаем с чистого состояния.
     */
    lu_init(lu);

    lu->n = n;

    lu->data = malloc(element_count * sizeof(double));
    lu->pivots = malloc(n * sizeof(size_t));

    if (lu->data == NULL || lu->pivots == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate LU factorization.\n");

        lu_free(lu);
        return -1;
    }

    /*
     * Исходную матрицу не разрушаем.
     */
    memcpy(
        lu->data,
        input,
        element_count * sizeof(double)
    );

    /*
     * Находим общий масштаб матрицы.
     */
    lu->matrix_scale = matrix_max_abs(
        input,
        n,
        options->parallel_threads
    );

    if (!isfinite(lu->matrix_scale) || lu->matrix_scale <= 0.0)
    {
        fprintf(
            stderr,
            "ERROR: LU input matrix is zero or contains invalid values.\n"
        );

        lu_free(lu);
        return -2;
    }

    double pivot_limit =
    options->pivot_tolerance * lu->matrix_scale;

    lu->swap_count = 0;


    /*
     * ========================================================
     * ОСНОВНОЙ ЦИКЛ LU
     * ========================================================
     *
     * Внешний цикл по k принципиально последовательный:
     *
     * шаг k+1 зависит от результата шага k.
     */
    for (size_t k = 0; k < n; ++k)
    {
        /*
         * ----------------------------------------------------
         * 1. ВЫБОР ГЛАВНОГО ЭЛЕМЕНТА
         * ----------------------------------------------------
         *
         * Partial pivoting:
         *
         *     max |A_ik|, i >= k.
         *
         * Этот цикл O(n), тогда как обновление хвоста O(n^2),
         * поэтому отдельный parallel region здесь обычно
         * не окупается.
         */
        size_t pivot_row = k;

        double pivot_abs =
        fabs(lu->data[k * n + k]);

        for (size_t i = k + 1; i < n; ++i)
        {
            double candidate =
            fabs(lu->data[i * n + k]);

            if (candidate > pivot_abs)
            {
                pivot_abs = candidate;
                pivot_row = i;
            }
        }

        /*
         * Проверяем численную вырожденность.
         */
        if (!isfinite(pivot_abs) || pivot_abs <= pivot_limit)
        {
            fprintf(
                stderr,
                "ERROR: LU matrix is singular or nearly singular "
                "at column %zu, |pivot| = %.9e.\n",
                k,
                pivot_abs
            );

            lu->factorized = 0;
            return -2;
        }

        /*
         * Сохраняем перестановку.
         *
         * В lu_solve() те же перестановки будут последовательно
         * применены к правой части b.
         */
        lu->pivots[k] = pivot_row;


        /*
         * ----------------------------------------------------
         * 2. ПЕРЕСТАНОВКА СТРОК
         * ----------------------------------------------------
         */
        if (pivot_row != k)
        {
            swap_rows(
                lu->data,
                n,
                k,
                pivot_row
            );

            ++lu->swap_count;
        }


        /*
         * Последний шаг уже не имеет хвостовой подматрицы.
         */
        if (k + 1 >= n)
        {
            continue;
        }

        double pivot = lu->data[k * n + k];

        /*
         * Дополнительная защита после перестановки.
         */
        if (!isfinite(pivot) || fabs(pivot) <= pivot_limit)
        {
            fprintf(
                stderr,
                "ERROR: invalid LU pivot at step %zu.\n",
                k
            );

            lu->factorized = 0;
            return -2;
        }


        /*
         * ----------------------------------------------------
         * 3-4. ELIMINATION + SCHUR UPDATE
         * ----------------------------------------------------
         *
         * Для каждой строки:
         *
         *     L_ik = A_ik / pivot
         *
         * и затем:
         *
         *     A_ij -= L_ik * A_kj.
         *
         *
         * Каждая строка i независима от других строк i,
         * поэтому эта часть хорошо распараллеливается.
         */
        size_t remaining = n - k - 1;

        #pragma omp parallel for schedule(static) num_threads(options->parallel_threads) if(remaining >= options->parallel_threshold)
        for (long long i_signed = (long long)k + 1;
             i_signed < (long long)n;
        ++i_signed)
             {
                 size_t i = (size_t)i_signed;

                 /*
                  * Множитель нижней треугольной матрицы:
                  *
                  *     L_ik.
                  *
                  * Он записывается прямо в нижнюю часть lu->data.
                  */
                 double multiplier =
                 lu->data[i * n + k] / pivot;

                 lu->data[i * n + k] = multiplier;


                 /*
                  * Обновляем строку хвостовой подматрицы.
                  *
                  * Внутренний цикл также хорошо векторизуется.
                  */
                 #pragma omp simd
                 for (size_t j = k + 1; j < n; ++j)
                 {
                     lu->data[i * n + j] -=
                     multiplier * lu->data[k * n + j];
                 }
             }
    }


    /*
     * Финальная проверка диагонали U.
     */
    for (size_t i = 0; i < n; ++i)
    {
        double diagonal =
        lu->data[i * n + i];

        if (!isfinite(diagonal) ||
            fabs(diagonal) <= pivot_limit)
        {
            fprintf(
                stderr,
                "ERROR: invalid diagonal U[%zu,%zu] = %.9e.\n",
                i,
                i,
                diagonal
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
 * ПРИМЕНЕНИЕ ПЕРЕСТАНОВОК К ПРАВОЙ ЧАСТИ
 * ============================================================
 *
 * Выполняем:
 *
 *     b <- P b.
 *
 * Важно применять перестановки в том же порядке,
 * в котором они выполнялись при факторизации.
 */
static void apply_pivots(
    const LUFactorization *lu,
    double *x
)
{
    for (size_t k = 0; k < lu->n; ++k)
    {
        size_t pivot_row = lu->pivots[k];

        if (pivot_row != k)
        {
            double temp = x[k];
            x[k] = x[pivot_row];
            x[pivot_row] = temp;
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
 * Диагональ L равна единице:
 *
 *     L_ii = 1.
 *
 * Поэтому деление на диагональ не требуется.
 *
 * Массив x содержит Pb на входе и y на выходе.
 */
static void forward_substitution(
    const LUFactorization *lu,
    double *x
)
{
    size_t n = lu->n;

    for (size_t i = 0; i < n; ++i)
    {
        double correction = 0.0;

        #pragma omp simd reduction(+:correction)
        for (size_t j = 0; j < i; ++j)
        {
            correction +=
            lu->data[i * n + j] * x[j];
        }

        x[i] -= correction;
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
static int backward_substitution(
    const LUFactorization *lu,
    double *x
)
{
    size_t n = lu->n;

    for (size_t ii = n; ii-- > 0;)
    {
        size_t i = ii;

        double correction = 0.0;

        #pragma omp simd reduction(+:correction)
        for (size_t j = i + 1; j < n; ++j)
        {
            correction +=
            lu->data[i * n + j] * x[j];
        }

        double diagonal =
        lu->data[i * n + i];

        if (diagonal == 0.0 || !isfinite(diagonal))
        {
            return -1;
        }

        x[i] =
        (x[i] - correction) / diagonal;
    }

    return 0;
}


/*
 * ============================================================
 * РЕШЕНИЕ ОДНОЙ СИСТЕМЫ
 * ============================================================
 */
int lu_solve(
    const LUFactorization *lu,
    const double *b,
    double *x
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

    size_t n = lu->n;

    /*
     * Если x и b разные массивы, копируем RHS.
     *
     * Если x == b, пользователь явно разрешил решение
     * "на месте", и копирование не требуется.
     */
    if (x != b)
    {
        memcpy(
            x,
            b,
            n * sizeof(double)
        );
    }

    /*
     * PA = LU
     *
     * A x = b
     *
     * => L U x = P b.
     */
    apply_pivots(lu, x);

    forward_substitution(lu, x);

    if (backward_substitution(lu, x) != 0)
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * РЕШЕНИЕ НЕСКОЛЬКИХ ПРАВЫХ ЧАСТЕЙ
 * ============================================================
 *
 * Для одной LU-факторизации часто потребуется решить систему
 * много раз.
 *
 * Например:
 *
 * - несколько направлений падающей волны;
 * - несколько поляризаций;
 * - несколько RHS для вспомогательных задач.
 *
 * Разные RHS полностью независимы.
 */
int lu_solve_multiple(
    const LUFactorization *lu,
    const double *B,
    double *X,
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
             size_t rhs = (size_t)rhs_signed;

             const double *b =
             &B[rhs * lu->n];

             double *x =
             &X[rhs * lu->n];

             int status = lu_solve(
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
 * ОПРЕДЕЛИТЕЛЬ
 * ============================================================
 */
double lu_determinant(const LUFactorization *lu)
{
    if (lu == NULL ||
        lu->data == NULL ||
        !lu->factorized)
    {
        return NAN;
    }

    /*
     * Каждая перестановка строк меняет знак det.
     */
    double determinant =
    (lu->swap_count % 2 == 0) ? 1.0 : -1.0;

    for (size_t i = 0; i < lu->n; ++i)
    {
        determinant *=
        lu->data[i * lu->n + i];
    }

    return determinant;
}


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void lu_print_info(const LUFactorization *lu)
{
    if (lu == NULL)
    {
        return;
    }

    printf("\n");
    printf("========================================\n");
    printf("LU FACTORIZATION\n");
    printf("========================================\n");
    printf("Size         : %zu x %zu\n", lu->n, lu->n);
    printf("Swap count   : %zu\n", lu->swap_count);
    printf("Matrix scale : %.9e\n", lu->matrix_scale);
    printf(
        "Factorized   : %s\n",
        lu->factorized ? "yes" : "no"
    );

    if (lu->data != NULL && lu->n > 0)
    {
        double memory_mib =
        ((double)lu->n
        * (double)lu->n
        * sizeof(double))
        / (1024.0 * 1024.0);

        printf("LU memory    : %.3f MiB\n", memory_mib);
    }

    printf("========================================\n");
}
