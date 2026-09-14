#include "harmonic.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void harmonic_system_init(
    HarmonicSystem *system
)
{
    if (system == NULL)
    {
        return;
    }

    memset(system, 0, sizeof(*system));

    complex_matrix_init(&system->matrix);
    complex_lu_init(&system->lu);
}


void harmonic_system_free(
    HarmonicSystem *system
)
{
    if (system == NULL)
    {
        return;
    }

    complex_matrix_free(&system->matrix);
    complex_lu_free(&system->lu);

    memset(system, 0, sizeof(*system));
}


void harmonic_solution_init(
    HarmonicSolution *solution
)
{
    if (solution == NULL)
    {
        return;
    }

    memset(solution, 0, sizeof(*solution));
}


void harmonic_solution_free(
    HarmonicSolution *solution
)
{
    if (solution == NULL)
    {
        return;
    }

    free(solution->edge_current);
    free(solution->node_voltage);

    harmonic_solution_init(solution);
}


/*
 * ============================================================
 * ПРОВЕРКА РАЗМЕРОВ
 * ============================================================
 */
static int harmonic_check_dimensions(
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const Matrix *PAT
)
{
    if (L == NULL ||
        R == NULL ||
        A == NULL ||
        PAT == NULL)
    {
        return -1;
    }

    if (L->data == NULL ||
        R->diagonal == NULL ||
        A->data == NULL ||
        PAT->data == NULL)
    {
        return -1;
    }

    size_t ne = L->n;
    size_t nv = A->cols;

    if (ne == 0 || nv == 0)
    {
        return -1;
    }

    if (A->rows != ne)
    {
        fprintf(
            stderr,
            "ERROR: harmonic system: A rows must equal Ne.\n"
        );

        return -1;
    }

    if (R->n != ne)
    {
        fprintf(
            stderr,
            "ERROR: harmonic system: R size must equal Ne.\n"
        );

        return -1;
    }

    if (PAT->rows != nv ||
        PAT->cols != ne)
    {
        fprintf(
            stderr,
            "ERROR: harmonic system: PAT must have size Nv x Ne.\n"
        );

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ПОСТРОЕНИЕ ПОЛНОЙ ГАРМОНИЧЕСКОЙ PEEC-МАТРИЦЫ
 * ============================================================
 *
 * Используем:
 *
 *     exp(j*w*t).
 *
 *
 * Система:
 *
 *     (R+j*w*L)I + A V = -U_e
 *
 *     -P A^T I + j*w V = P I_s.
 *
 *
 * Матрица:
 *
 *     [ R+j*w*L      A     ]
 *     [                    ]
 *     [ -P*A^T      j*w*Iv ].
 */
int harmonic_build_system(
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const Matrix *PAT,
    double frequency,
    int parallel_threads,
    HarmonicSystem *system
)
{
    if (system == NULL ||
        parallel_threads < 1 ||
        !isfinite(frequency) ||
        frequency <= 0.0)
    {
        return -1;
    }

    if (harmonic_check_dimensions(
            L,
            R,
            A,
            PAT
        ) != 0)
    {
        return -1;
    }

    harmonic_system_free(system);
    harmonic_system_init(system);

    size_t ne = L->n;
    size_t nv = A->cols;
    size_t n = ne + nv;

    double omega =
        2.0 * M_PI * frequency;

    system->n_edges = ne;
    system->n_nodes = nv;
    system->system_size = n;
    system->frequency = frequency;
    system->omega = omega;

    if (complex_matrix_alloc(
            &system->matrix,
            n,
            n
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate harmonic PEEC matrix.\n"
        );

        harmonic_system_free(system);

        return -1;
    }


    /*
     * ========================================================
     * ВЕРХНИЙ ЛЕВЫЙ БЛОК
     * ========================================================
     *
     *     R + j*w*L.
     *
     * Размер:
     *
     *     Ne x Ne.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)ne;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        for (size_t j = 0; j < ne; ++j)
        {
            Complex value;

            value.re =
                (i == j)
                ? R->diagonal[i]
                : 0.0;

            value.im =
                omega
                * L->data[i * ne + j];

            system->matrix.data[i * n + j] =
                value;
        }
    }


    /*
     * ========================================================
     * ВЕРХНИЙ ПРАВЫЙ БЛОК
     * ========================================================
     *
     *     A.
     *
     * Размер:
     *
     *     Ne x Nv.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)ne;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        for (size_t j = 0; j < nv; ++j)
        {
            size_t row = i;
            size_t col = ne + j;

            system->matrix.data[row * n + col].re =
                A->data[i * nv + j];

            system->matrix.data[row * n + col].im =
                0.0;
        }
    }


    /*
     * ========================================================
     * НИЖНИЙ ЛЕВЫЙ БЛОК
     * ========================================================
     *
     *     -P A^T.
     *
     * Размер:
     *
     *     Nv x Ne.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)nv;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        for (size_t j = 0; j < ne; ++j)
        {
            size_t row = ne + i;
            size_t col = j;

            system->matrix.data[row * n + col].re =
                -PAT->data[i * ne + j];

            system->matrix.data[row * n + col].im =
                0.0;
        }
    }


    /*
     * ========================================================
     * НИЖНИЙ ПРАВЫЙ БЛОК
     * ========================================================
     *
     *     j*w*Iv.
     *
     * Поскольку complex_matrix_alloc() должна
     * инициализировать память нулями, достаточно
     * записать только диагональ.
     *
     * Если твоя complex_matrix_alloc() НЕ обнуляет
     * память, перед заполнением этой матрицы нужно
     * вызвать complex_matrix_zero().
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)nv;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        size_t row = ne + i;
        size_t col = ne + i;

        system->matrix.data[row * n + col].re =0.0;

        system->matrix.data[row * n + col].im =omega;
    }


    system->matrix_built = 1;
    system->factorized = 0;

    return 0;
}


/*
 * ============================================================
 * LU-ФАКТОРИЗАЦИЯ
 * ============================================================
 */
int harmonic_factorize(
    HarmonicSystem *system,
    int parallel_threads
)
{
    if (system == NULL ||
        !system->matrix_built ||
        system->matrix.data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    complex_lu_free(&system->lu);
    complex_lu_init(&system->lu);

    ComplexLUOptions options;

    complex_lu_options_default(
        &options,
        parallel_threads
    );

    if (complex_lu_factorize(
            system->matrix.data,
            system->system_size,
            &options,
            &system->lu
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic PEEC LU factorization failed.\n"
        );

        system->factorized = 0;

        return -1;
    }

    system->factorized = 1;

    return 0;
}


/*
 * ============================================================
 * ПРАВАЯ ЧАСТЬ
 * ============================================================
 *
 *     rhs =
 *
 *     [ -U_e  ]
 *     [ P I_s ].
 */
int harmonic_build_rhs(
    const PotentialMatrix *P,
    const HarmonicExcitation *excitation,
    Complex *rhs,
    int parallel_threads
)
{
    if (P == NULL ||
        excitation == NULL ||
        rhs == NULL ||
        P->data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    size_t ne = excitation->n_edges;
    size_t nv = excitation->n_nodes;

    if (P->n != nv)
    {
        fprintf(
            stderr,
            "ERROR: harmonic RHS: P size does not match Nv.\n"
        );

        return -1;
    }


    /*
     * ========================================================
     * ВЕРХНЯЯ ЧАСТЬ
     * ========================================================
     *
     *     -U_e.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)ne;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        if (excitation->edge_voltage != NULL)
        {
            rhs[i].re =
                -excitation->edge_voltage[i].re;

            rhs[i].im =
                -excitation->edge_voltage[i].im;
        }
        else
        {
            rhs[i].re = 0.0;
            rhs[i].im = 0.0;
        }
    }


    /*
     * ========================================================
     * НИЖНЯЯ ЧАСТЬ
     * ========================================================
     *
     *     P I_s.
     */
    if (excitation->node_current == NULL)
    {
        #pragma omp parallel for schedule(static) num_threads(parallel_threads)
        for (long long i_signed = 0;
             i_signed < (long long)nv;
             ++i_signed)
        {
            size_t i = (size_t)i_signed;

            rhs[ne + i].re = 0.0;
            rhs[ne + i].im = 0.0;
        }

        return 0;
    }


    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)nv;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        double sum_re = 0.0;
        double sum_im = 0.0;

        const double *p_row =
            &P->data[i * nv];

        #pragma omp simd reduction(+:sum_re,sum_im)
        for (size_t j = 0; j < nv; ++j)
        {
            double p = p_row[j];

            sum_re +=
                p
                * excitation->node_current[j].re;

            sum_im +=
                p
                * excitation->node_current[j].im;
        }

        rhs[ne + i].re = sum_re;
        rhs[ne + i].im = sum_im;
    }

    return 0;
}


/*
 * ============================================================
 * РЕШЕНИЕ
 * ============================================================
 */
int harmonic_solve(
    const HarmonicSystem *system,
    const PotentialMatrix *P,
    const HarmonicExcitation *excitation,
    HarmonicSolution *solution,
    int parallel_threads
)
{
    if (system == NULL ||
        P == NULL ||
        excitation == NULL ||
        solution == NULL ||
        !system->factorized ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (excitation->n_edges != system->n_edges ||
        excitation->n_nodes != system->n_nodes)
    {
        fprintf(
            stderr,
            "ERROR: harmonic excitation dimensions do not match system.\n"
        );

        return -1;
    }

    size_t ne = system->n_edges;
    size_t nv = system->n_nodes;
    size_t n = system->system_size;

    Complex *rhs =
        malloc(n * sizeof(Complex));

    Complex *x =
        malloc(n * sizeof(Complex));

    if (rhs == NULL ||
        x == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate harmonic solution vectors.\n"
        );

        free(rhs);
        free(x);

        return -1;
    }


    if (harmonic_build_rhs(
            P,
            excitation,
            rhs,
            parallel_threads
        ) != 0)
    {
        free(rhs);
        free(x);

        return -1;
    }


    if (complex_lu_solve(
            &system->lu,
            rhs,
            x
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic PEEC solve failed.\n"
        );

        free(rhs);
        free(x);

        return -1;
    }


    harmonic_solution_free(solution);

    solution->n_edges = ne;
    solution->n_nodes = nv;

    solution->edge_current =
        malloc(ne * sizeof(Complex));

    solution->node_voltage =
        malloc(nv * sizeof(Complex));

    if (solution->edge_current == NULL ||
        solution->node_voltage == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate harmonic output.\n"
        );

        harmonic_solution_free(solution);

        free(rhs);
        free(x);

        return -1;
    }


    memcpy(
        solution->edge_current,
        x,
        ne * sizeof(Complex)
    );

    memcpy(
        solution->node_voltage,
        &x[ne],
        nv * sizeof(Complex)
    );


    free(rhs);
    free(x);

    return 0;
}


/*
 * ============================================================
 * НЕВЯЗКА ПОЛНОЙ PEEC-СИСТЕМЫ
 * ============================================================
 *
 * Проверяем:
 *
 *     Z x = b.
 *
 *
 * Для каждой строки:
 *
 *     r_i = (Zx-b)_i.
 *
 *
 * Дополнительно считаем:
 *
 *                       |r_i|
 *     eta_i = ----------------------------
 *             sum_j |Z_ij||x_j| + |b_i|.
 *
 *
 * Итоговая backward error:
 *
 *     eta = max_i eta_i.
 */
int harmonic_compute_residual(
    const HarmonicSystem *system,
    const PotentialMatrix *P,
    const HarmonicExcitation *excitation,
    const HarmonicSolution *solution,
    int parallel_threads,
    HarmonicResidual *residual
)
{
    if (system == NULL ||
        P == NULL ||
        excitation == NULL ||
        solution == NULL ||
        residual == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (!system->matrix_built ||
        system->matrix.data == NULL)
    {
        fprintf(
            stderr,
            "ERROR: harmonic residual requires a built system matrix.\n"
        );

        return -1;
    }

    if (excitation->n_edges != system->n_edges ||
        excitation->n_nodes != system->n_nodes)
    {
        fprintf(
            stderr,
            "ERROR: harmonic residual excitation dimensions mismatch.\n"
        );

        return -1;
    }

    if (solution->n_edges != system->n_edges ||
        solution->n_nodes != system->n_nodes ||
        solution->edge_current == NULL ||
        solution->node_voltage == NULL)
    {
        fprintf(
            stderr,
            "ERROR: harmonic residual solution dimensions mismatch.\n"
        );

        return -1;
    }

    size_t ne = system->n_edges;
    size_t nv = system->n_nodes;
    size_t n = system->system_size;

    Complex *x =
        malloc(n * sizeof(Complex));

    Complex *rhs =
        malloc(n * sizeof(Complex));

    Complex *y =
        malloc(n * sizeof(Complex));

    if (x == NULL ||
        rhs == NULL ||
        y == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate harmonic residual vectors.\n"
        );

        free(x);
        free(rhs);
        free(y);

        return -1;
    }


    /*
     * ========================================================
     * x = [ I ; V ]
     * ========================================================
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)ne;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        x[i] =
            solution->edge_current[i];
    }


    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)nv;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        x[ne + i] =
            solution->node_voltage[i];
    }


    /*
     * ========================================================
     * b = [ -U_e ; P I_s ]
     * ========================================================
     */
    if (harmonic_build_rhs(
            P,
            excitation,
            rhs,
            parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot build RHS for harmonic residual.\n"
        );

        free(x);
        free(rhs);
        free(y);

        return -1;
    }


    /*
     * ========================================================
     * y = Z x
     * ========================================================
     */
    if (complex_matrix_vector_multiply(
            &system->matrix,
            x,
            y,
            parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic residual matrix-vector multiply failed.\n"
        );

        free(x);
        free(rhs);
        free(y);

        return -1;
    }


    double absolute_residual = 0.0;
    double rhs_inf = 0.0;

    double edge_absolute_residual = 0.0;
    double node_absolute_residual = 0.0;

    double backward_error = 0.0;
    double edge_backward_error = 0.0;
    double node_backward_error = 0.0;


    /*
     * ========================================================
     * ПОСТРОЧНАЯ ДИАГНОСТИКА
     * ========================================================
     */
    #pragma omp parallel for \
        reduction(max:absolute_residual,rhs_inf,edge_absolute_residual,node_absolute_residual,backward_error,edge_backward_error,node_backward_error) \
        schedule(static) num_threads(parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)n;
         ++i_signed)
    {
        size_t i = (size_t)i_signed;

        double residual_re =
            y[i].re - rhs[i].re;

        double residual_im =
            y[i].im - rhs[i].im;

        double residual_abs =
            hypot(
                residual_re,
                residual_im
            );

        double rhs_abs =
            hypot(
                rhs[i].re,
                rhs[i].im
            );


        /*
         * Нормировочный знаменатель:
         *
         *     sum_j |Z_ij||x_j| + |b_i|.
         */
        double denominator =
            rhs_abs;

        const Complex *row =
            &system->matrix.data[i * n];


        for (size_t j = 0; j < n; ++j)
        {
            double zij_abs =
                hypot(
                    row[j].re,
                    row[j].im
                );

            double xj_abs =
                hypot(
                    x[j].re,
                    x[j].im
                );

            denominator +=
                zij_abs * xj_abs;
        }


        double eta = 0.0;

        if (denominator > 0.0)
        {
            eta =
                residual_abs
                / denominator;
        }
        else
        {
            eta =
                residual_abs;
        }


        if (residual_abs > absolute_residual)
        {
            absolute_residual =
                residual_abs;
        }

        if (rhs_abs > rhs_inf)
        {
            rhs_inf =
                rhs_abs;
        }

        if (eta > backward_error)
        {
            backward_error =
                eta;
        }


        /*
         * ====================================================
         * ВЕРХНИЙ PEEC-БЛОК
         * ====================================================
         *
         *     (R+j*w*L)I + AV = -U_e.
         */
        if (i < ne)
        {
            if (residual_abs >
                edge_absolute_residual)
            {
                edge_absolute_residual =
                    residual_abs;
            }

            if (eta >
                edge_backward_error)
            {
                edge_backward_error =
                    eta;
            }
        }


        /*
         * ====================================================
         * НИЖНИЙ PEEC-БЛОК
         * ====================================================
         *
         *     -P A^T I + j*w V = P I_s.
         */
        else
        {
            if (residual_abs >
                node_absolute_residual)
            {
                node_absolute_residual =
                    residual_abs;
            }

            if (eta >
                node_backward_error)
            {
                node_backward_error =
                    eta;
            }
        }
    }


    residual->absolute_residual =
        absolute_residual;

    residual->edge_absolute_residual =
        edge_absolute_residual;

    residual->node_absolute_residual =
        node_absolute_residual;

    residual->backward_error =
        backward_error;

    residual->edge_backward_error =
        edge_backward_error;

    residual->node_backward_error =
        node_backward_error;


    if (rhs_inf > 0.0)
    {
        residual->rhs_relative_residual =
            absolute_residual
            / rhs_inf;
    }
    else
    {
        residual->rhs_relative_residual =
            absolute_residual;
    }


    free(x);
    free(rhs);
    free(y);

    return 0;
}


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void harmonic_print_info(
    const HarmonicSystem *system
)
{
    if (system == NULL)
    {
        return;
    }

    printf("\n");
    printf("========================================\n");
    printf("HARMONIC PEEC SYSTEM\n");
    printf("========================================\n");

    printf(
        "Frequency          : %.9e Hz\n",
        system->frequency
    );

    printf(
        "Angular frequency  : %.9e rad/s\n",
        system->omega
    );

    printf(
        "Edge unknowns      : %zu\n",
        system->n_edges
    );

    printf(
        "Node unknowns      : %zu\n",
        system->n_nodes
    );

    printf(
        "Full system size   : %zu x %zu\n",
        system->system_size,
        system->system_size
    );

    if (system->matrix.data != NULL)
    {
        double memory_mib =
            (double)system->system_size
            * (double)system->system_size
            * sizeof(Complex)
            / (1024.0 * 1024.0);

        printf(
            "System memory      : %.3f MiB\n",
            memory_mib
        );
    }

    printf(
        "Matrix built       : %s\n",
        system->matrix_built
        ? "yes"
        : "no"
    );

    printf(
        "LU factorized      : %s\n",
        system->factorized
        ? "yes"
        : "no"
    );

    printf("========================================\n");
}