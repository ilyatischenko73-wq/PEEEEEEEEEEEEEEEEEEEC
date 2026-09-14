#include "config.h"
#include "mesh.h"
#include "topology.h"
#include "dual_mesh.h"
#include "incidence.h"
#include "inductance.h"
#include "potential.h"
#include "resistance.h"
#include "matrix.h"
#include "harmonic.h"
#include "incident_field.h"
#include "scattering.h"

#ifdef TEST
#include "lu.h"
#include "complex_matrix.h"
#include "complex_lu.h"
#include "transient.h"
#endif

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


/*
 * ============================================================
 * ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
 * ============================================================
 */

static const char *physics_name(PeecPhysicsMode mode)
{
    switch (mode)
    {
        case PEEC_PHYSICS_QUASISTATIC:
            return "quasistatic";

        case PEEC_PHYSICS_RETARDED:
            return "retarded";

        default:
            return "unknown";
    }
}


static const char *task_name(PeecTask task)
{
    switch (task)
    {
        case PEEC_TASK_MESH_INFO:
            return "mesh-info";

        case PEEC_TASK_SCATTERING:
            return "scattering";

        case PEEC_TASK_LIGHTNING:
            return "lightning";

        case PEEC_TASK_NONE:
        default:
            return "none";
    }
}


/*
 * ============================================================
 * МАКСИМАЛЬНЫЙ МОДУЛЬ КОМПЛЕКСНОГО ВЕКТОРА
 * ============================================================
 */
static double complex_vector_max_abs(
    const Complex *x,
    size_t n
)
{
    if (x == NULL || n == 0)
    {
        return 0.0;
    }

    double max_value = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        double value = hypot(x[i].re, x[i].im);

        if (value > max_value)
        {
            max_value = value;
        }
    }

    return max_value;
}


#ifdef TEST

/*
 * ============================================================
 * TEST 1. ВЕЩЕСТВЕННЫЙ LU
 * ============================================================
 */
static int test_lu(int parallel_threads)
{
    const size_t n = 3;
    const double tolerance = 1.0e-12;

    double A[9] = {
        4.0, 2.0, 1.0,
        2.0, 5.0, 2.0,
        1.0, 2.0, 4.0
    };

    double x_true[3] = {
        1.0,
        2.0,
        3.0
    };

    double b[3] = {
        11.0,
        18.0,
        17.0
    };

    double x[3] = {
        0.0,
        0.0,
        0.0
    };

    LUOptions options;
    LUFactorization lu;

    lu_options_default(&options, parallel_threads);
    lu_init(&lu);

    printf("\n");
    printf("============================================================\n");
    printf("REAL LU SELF-TEST\n");
    printf("============================================================\n");

    if (lu_factorize(A, n, &options, &lu) != 0)
    {
        fprintf(stderr, "ERROR: real LU factorization test failed.\n");
        lu_free(&lu);
        return -1;
    }

    if (lu_solve(&lu, b, x) != 0)
    {
        fprintf(stderr, "ERROR: real LU solve test failed.\n");
        lu_free(&lu);
        return -1;
    }

    double max_solution_error = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        double error = fabs(x[i] - x_true[i]);

        if (error > max_solution_error)
        {
            max_solution_error = error;
        }
    }

    double max_residual = 0.0;
    double b_norm = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        double ax = 0.0;

        for (size_t j = 0; j < n; ++j)
        {
            ax += A[i * n + j] * x[j];
        }

        double residual = fabs(ax - b[i]);

        if (residual > max_residual)
        {
            max_residual = residual;
        }

        double value = fabs(b[i]);

        if (value > b_norm)
        {
            b_norm = value;
        }
    }

    double relative_residual = 0.0;

    if (b_norm > 0.0)
    {
        relative_residual = max_residual / b_norm;
    }

    for (size_t i = 0; i < n; ++i)
    {
        printf(
            "x[%zu] = %.16e   expected %.16e\n",
            i,
            x[i],
            x_true[i]
        );
    }

    printf("\n");
    printf("Max solution error : %.3e\n", max_solution_error);
    printf("Max residual       : %.3e\n", max_residual);
    printf("Relative residual  : %.3e\n", relative_residual);

    if (max_solution_error > tolerance ||
        relative_residual > tolerance)
    {
        fprintf(stderr, "ERROR: real LU self-test failed.\n");
        lu_free(&lu);
        return -1;
    }

    printf("Real LU test        : PASSED\n");
    printf("============================================================\n");

    lu_free(&lu);

    return 0;
}


/*
 * ============================================================
 * TEST 2. КОМПЛЕКСНАЯ АРИФМЕТИКА
 * ============================================================
 */
static int test_complex_arithmetic(void)
{
    const double tolerance = 1.0e-14;

    Complex a = complex_make(1.0, 2.0);
    Complex b = complex_make(3.0, -4.0);

    Complex sum = complex_add(a, b);
    Complex difference = complex_sub(a, b);
    Complex product = complex_mul(a, b);
    Complex quotient = complex_div(a, b);
    Complex conjugate = complex_conj(a);

    double abs_a = complex_abs(a);

    printf("\n");
    printf("============================================================\n");
    printf("COMPLEX ARITHMETIC SELF-TEST\n");
    printf("============================================================\n");

    if (fabs(sum.re - 4.0) > tolerance ||
        fabs(sum.im + 2.0) > tolerance)
    {
        fprintf(stderr, "ERROR: complex addition test failed.\n");
        return -1;
    }

    if (fabs(difference.re + 2.0) > tolerance ||
        fabs(difference.im - 6.0) > tolerance)
    {
        fprintf(stderr, "ERROR: complex subtraction test failed.\n");
        return -1;
    }

    if (fabs(product.re - 11.0) > tolerance ||
        fabs(product.im - 2.0) > tolerance)
    {
        fprintf(stderr, "ERROR: complex multiplication test failed.\n");
        return -1;
    }

    if (fabs(quotient.re + 0.2) > tolerance ||
        fabs(quotient.im - 0.4) > tolerance)
    {
        fprintf(stderr, "ERROR: complex division test failed.\n");
        return -1;
    }

    if (fabs(conjugate.re - 1.0) > tolerance ||
        fabs(conjugate.im + 2.0) > tolerance)
    {
        fprintf(stderr, "ERROR: complex conjugation test failed.\n");
        return -1;
    }

    if (fabs(abs_a - sqrt(5.0)) > tolerance)
    {
        fprintf(stderr, "ERROR: complex absolute-value test failed.\n");
        return -1;
    }

    printf("a + b   = %.16e %+.16ej\n", sum.re, sum.im);
    printf("a - b   = %.16e %+.16ej\n", difference.re, difference.im);
    printf("a * b   = %.16e %+.16ej\n", product.re, product.im);
    printf("a / b   = %.16e %+.16ej\n", quotient.re, quotient.im);
    printf("conj(a) = %.16e %+.16ej\n", conjugate.re, conjugate.im);
    printf("|a|     = %.16e\n", abs_a);
    printf("Complex arithmetic test : PASSED\n");
    printf("============================================================\n");

    return 0;
}


/*
 * ============================================================
 * TEST 3. COMPLEX MATRIX
 * ============================================================
 */
static int test_complex_matrix(int parallel_threads)
{
    const double tolerance = 1.0e-13;

    ComplexMatrix A;
    ComplexMatrix AT;
    ComplexMatrix AH;

    complex_matrix_init(&A);
    complex_matrix_init(&AT);
    complex_matrix_init(&AH);

    if (complex_matrix_alloc(&A, 2, 2) != 0)
    {
        return -1;
    }

    A.data[0] = complex_make(1.0, 1.0);
    A.data[1] = complex_make(2.0, 0.0);
    A.data[2] = complex_make(3.0, -1.0);
    A.data[3] = complex_make(-1.0, 2.0);

    Complex x[2];

    x[0] = complex_make(1.0, -1.0);
    x[1] = complex_make(2.0, 1.0);

    Complex y[2];

    y[0] = complex_make(0.0, 0.0);
    y[1] = complex_make(0.0, 0.0);

    printf("\n");
    printf("============================================================\n");
    printf("COMPLEX MATRIX SELF-TEST\n");
    printf("============================================================\n");

    if (complex_matrix_vector_multiply(
            &A,
            x,
            y,
            parallel_threads
        ) != 0)
    {
        complex_matrix_free(&A);
        complex_matrix_free(&AT);
        complex_matrix_free(&AH);

        return -1;
    }

    Complex expected_y[2];

    expected_y[0] = complex_make(6.0, 2.0);
    expected_y[1] = complex_make(-2.0, -1.0);

    double max_error = 0.0;

    for (size_t i = 0; i < 2; ++i)
    {
        double error_re = y[i].re - expected_y[i].re;
        double error_im = y[i].im - expected_y[i].im;
        double error = hypot(error_re, error_im);

        if (error > max_error)
        {
            max_error = error;
        }
    }

    if (max_error > tolerance)
    {
        fprintf(stderr, "ERROR: complex matrix self-test failed.\n");

        complex_matrix_free(&A);
        complex_matrix_free(&AT);
        complex_matrix_free(&AH);

        return -1;
    }

    if (complex_matrix_transpose(
            &A,
            &AT,
            parallel_threads
        ) != 0)
    {
        fprintf(stderr, "ERROR: complex transpose self-test failed.\n");

        complex_matrix_free(&A);
        complex_matrix_free(&AT);
        complex_matrix_free(&AH);

        return -1;
    }

    Complex at01 = complex_matrix_get(&AT, 0, 1);

    if (fabs(at01.re - 3.0) > tolerance ||
        fabs(at01.im + 1.0) > tolerance)
    {
        fprintf(stderr, "ERROR: complex transpose value is incorrect.\n");

        complex_matrix_free(&A);
        complex_matrix_free(&AT);
        complex_matrix_free(&AH);

        return -1;
    }

    if (complex_matrix_conjugate_transpose(
            &A,
            &AH,
            parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: conjugate transpose self-test failed.\n"
        );

        complex_matrix_free(&A);
        complex_matrix_free(&AT);
        complex_matrix_free(&AH);

        return -1;
    }

    Complex ah01 = complex_matrix_get(&AH, 0, 1);

    if (fabs(ah01.re - 3.0) > tolerance ||
        fabs(ah01.im - 1.0) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: conjugate transpose value is incorrect.\n"
        );

        complex_matrix_free(&A);
        complex_matrix_free(&AT);
        complex_matrix_free(&AH);

        return -1;
    }

    printf(
        "y[0] = %.16e %+.16ej   expected %.16e %+.16ej\n",
        y[0].re,
        y[0].im,
        expected_y[0].re,
        expected_y[0].im
    );

    printf(
        "y[1] = %.16e %+.16ej   expected %.16e %+.16ej\n",
        y[1].re,
        y[1].im,
        expected_y[1].re,
        expected_y[1].im
    );

    printf("Max matrix-vector error : %.3e\n", max_error);
    printf("Complex matrix test     : PASSED\n");
    printf("============================================================\n");

    complex_matrix_free(&A);
    complex_matrix_free(&AT);
    complex_matrix_free(&AH);

    return 0;
}


/*
 * ============================================================
 * TEST 4. COMPLEX LU
 * ============================================================
 */
static int test_complex_lu(int parallel_threads)
{
    const size_t n = 3;
    const double tolerance = 1.0e-12;

    ComplexMatrix A;

    complex_matrix_init(&A);

    if (complex_matrix_alloc(&A, n, n) != 0)
    {
        return -1;
    }

    A.data[0] = complex_make(0.0, 0.0);
    A.data[1] = complex_make(2.0, 1.0);
    A.data[2] = complex_make(1.0, -1.0);

    A.data[3] = complex_make(3.0, 1.0);
    A.data[4] = complex_make(1.0, 0.0);
    A.data[5] = complex_make(0.0, -2.0);

    A.data[6] = complex_make(1.0, -1.0);
    A.data[7] = complex_make(4.0, 2.0);
    A.data[8] = complex_make(5.0, 0.0);

    Complex x_true[3];

    x_true[0] = complex_make(1.0, 2.0);
    x_true[1] = complex_make(-1.0, 1.0);
    x_true[2] = complex_make(2.0, -1.0);

    Complex b[3];

    if (complex_matrix_vector_multiply(
            &A,
            x_true,
            b,
            parallel_threads
        ) != 0)
    {
        complex_matrix_free(&A);
        return -1;
    }

    Complex x[3];

    for (size_t i = 0; i < n; ++i)
    {
        x[i] = complex_make(0.0, 0.0);
    }

    ComplexLUOptions options;
    ComplexLUFactorization lu;

    complex_lu_options_default(&options, parallel_threads);
    complex_lu_init(&lu);

    printf("\n");
    printf("============================================================\n");
    printf("COMPLEX LU SELF-TEST\n");
    printf("============================================================\n");

    if (complex_lu_factorize(
            A.data,
            n,
            &options,
            &lu
        ) != 0)
    {
        complex_lu_free(&lu);
        complex_matrix_free(&A);

        return -1;
    }

    if (lu.swap_count == 0)
    {
        fprintf(
            stderr,
            "ERROR: complex LU pivoting self-test failed.\n"
        );

        complex_lu_free(&lu);
        complex_matrix_free(&A);

        return -1;
    }

    if (complex_lu_solve(
            &lu,
            b,
            x
        ) != 0)
    {
        complex_lu_free(&lu);
        complex_matrix_free(&A);

        return -1;
    }

    double max_solution_error = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        double error_re = x[i].re - x_true[i].re;
        double error_im = x[i].im - x_true[i].im;
        double error = hypot(error_re, error_im);

        if (error > max_solution_error)
        {
            max_solution_error = error;
        }
    }

    Complex ax[3];

    if (complex_matrix_vector_multiply(
            &A,
            x,
            ax,
            parallel_threads
        ) != 0)
    {
        complex_lu_free(&lu);
        complex_matrix_free(&A);

        return -1;
    }

    double max_residual = 0.0;
    double b_norm = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        double residual_re = ax[i].re - b[i].re;
        double residual_im = ax[i].im - b[i].im;
        double residual = hypot(residual_re, residual_im);

        if (residual > max_residual)
        {
            max_residual = residual;
        }

        double b_abs = hypot(b[i].re, b[i].im);

        if (b_abs > b_norm)
        {
            b_norm = b_abs;
        }
    }

    double relative_residual = 0.0;

    if (b_norm > 0.0)
    {
        relative_residual = max_residual / b_norm;
    }

    for (size_t i = 0; i < n; ++i)
    {
        printf(
            "x[%zu] = %.16e %+.16ej   expected %.16e %+.16ej\n",
            i,
            x[i].re,
            x[i].im,
            x_true[i].re,
            x_true[i].im
        );
    }

    printf("\n");
    printf("Row swaps          : %zu\n", lu.swap_count);
    printf("Max solution error : %.3e\n", max_solution_error);
    printf("Max residual       : %.3e\n", max_residual);
    printf("Relative residual  : %.3e\n", relative_residual);

    if (max_solution_error > tolerance ||
        relative_residual > tolerance)
    {
        fprintf(stderr, "ERROR: complex LU self-test failed.\n");

        complex_lu_free(&lu);
        complex_matrix_free(&A);

        return -1;
    }

    printf("Complex LU test    : PASSED\n");
    printf("============================================================\n");

    complex_lu_free(&lu);
    complex_matrix_free(&A);

    return 0;
}


/*
 * ============================================================
 * TEST 5. СОПРОТИВЛЕНИЕ
 * ============================================================
 */
static int test_resistance(int parallel_threads)
{
    const double tolerance = 1.0e-12;

    Mesh mesh;
    DualMesh dual;

    mesh_init(&mesh);
    dual_mesh_init(&dual);

    mesh.n_edges = 2;

    mesh.edge_lengths = calloc(
        mesh.n_edges,
        sizeof(double)
    );

    if (mesh.edge_lengths == NULL)
    {
        return -1;
    }

    mesh.edge_lengths[0] = 2.0;
    mesh.edge_lengths[1] = 3.0;

    dual.n_edges = 2;

    dual.edge_region_areas = calloc(
        dual.n_edges,
        sizeof(double)
    );

    if (dual.edge_region_areas == NULL)
    {
        mesh_free(&mesh);
        return -1;
    }

    dual.edge_region_areas[0] = 1.0;
    dual.edge_region_areas[1] = 1.5;

    ResistanceOptions options;
    ResistanceMatrix R;

    resistance_options_default(&options, parallel_threads);
    resistance_init(&R);

    printf("\n");
    printf("============================================================\n");
    printf("RESISTANCE SELF-TEST\n");
    printf("============================================================\n");

    if (resistance_compute_matrix(
            &mesh,
            &dual,
            &options,
            &R
        ) != 0)
    {
        resistance_free(&R);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return -1;
    }

    for (size_t i = 0; i < R.n; ++i)
    {
        if (fabs(R.diagonal[i]) > tolerance)
        {
            resistance_free(&R);
            dual_mesh_free(&dual);
            mesh_free(&mesh);

            return -1;
        }
    }

    printf("PEC resistance test   : PASSED\n");

    resistance_free(&R);

    options.model = PEEC_CONDUCTOR_SHEET_RESISTANCE;
    options.sheet_resistance = 5.0;

    if (resistance_compute_matrix(
            &mesh,
            &dual,
            &options,
            &R
        ) != 0)
    {
        resistance_free(&R);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return -1;
    }

    if (fabs(R.diagonal[0] - 20.0) > tolerance ||
        fabs(R.diagonal[1] - 30.0) > tolerance)
    {
        resistance_free(&R);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return -1;
    }

    printf("Sheet resistance test : PASSED\n");
    printf("============================================================\n");

    resistance_free(&R);
    dual_mesh_free(&dual);
    mesh_free(&mesh);

    return 0;
}


/*
 * ============================================================
 * TEST 6. P A^T
 * ============================================================
 */
static int test_incidence_PAT(int parallel_threads)
{
    const double tolerance = 1.0e-14;

    Mesh mesh;
    PotentialMatrix P;
    Matrix PAT;

    mesh_init(&mesh);
    potential_init(&P);
    matrix_init(&PAT);

    mesh.n_nodes = 3;
    mesh.n_edges = 2;

    mesh.edges = malloc(
        2 * mesh.n_edges * sizeof(int)
    );

    if (mesh.edges == NULL)
    {
        return -1;
    }

    mesh.edges[0] = 0;
    mesh.edges[1] = 1;

    mesh.edges[2] = 1;
    mesh.edges[3] = 2;

    P.n = 3;

    P.data = malloc(
        P.n * P.n * sizeof(double)
    );

    if (P.data == NULL)
    {
        mesh_free(&mesh);
        return -1;
    }

    P.data[0] = 10.0;
    P.data[1] = 2.0;
    P.data[2] = 3.0;

    P.data[3] = 2.0;
    P.data[4] = 20.0;
    P.data[5] = 4.0;

    P.data[6] = 3.0;
    P.data[7] = 4.0;
    P.data[8] = 30.0;

    printf("\n");
    printf("============================================================\n");
    printf("P A^T SELF-TEST\n");
    printf("============================================================\n");

    if (incidence_build_PAT(
            &mesh,
            &P,
            &PAT,
            parallel_threads
        ) != 0)
    {
        matrix_free(&PAT);
        potential_free(&P);
        mesh_free(&mesh);

        return -1;
    }

    double expected[6] = {
        -8.0, 1.0,
        18.0, -16.0,
        1.0, 26.0
    };

    double max_error = 0.0;

    for (size_t i = 0; i < 6; ++i)
    {
        double error = fabs(
            PAT.data[i] - expected[i]
        );

        if (error > max_error)
        {
            max_error = error;
        }
    }

    printf("Max P A^T error : %.3e\n", max_error);

    if (max_error > tolerance)
    {
        matrix_free(&PAT);
        potential_free(&P);
        mesh_free(&mesh);

        return -1;
    }

    printf("P A^T test       : PASSED\n");
    printf("============================================================\n");

    matrix_free(&PAT);
    potential_free(&P);
    mesh_free(&mesh);

    return 0;
}



/*
 * ============================================================
 * TEST 7. TRANSIENT PEEC
 * ============================================================
 *
 * Используем минимальную PEEC-модель:
 *
 *     Ne = 1
 *     Nv = 2.
 *
 *
 * Единственное ребро:
 *
 *     node 0 --------> node 1
 *
 *
 * Поэтому:
 *
 *     A = [ -1  +1 ].
 *
 *
 * Выбираем:
 *
 *     L = 2 H
 *     R = 0
 *
 *     P = [ 1  0 ]
 *         [ 0  1 ].
 *
 *
 * Тогда:
 *
 *     P A^T =
 *
 *     [ -1 ]
 *     [ +1 ].
 *
 *
 * Для Stage 1 без источников:
 *
 *     2 dI/dt - V0 + V1 = 0
 *
 *     dV0/dt + I = 0
 *
 *     dV1/dt - I = 0.
 *
 *
 * При:
 *
 *     I(0)  = 1
 *     V0(0) = 0
 *     V1(0) = 0
 *
 * точное решение:
 *
 *     I(t)  = cos(t)
 *
 *     V0(t) = -sin(t)
 *
 *     V1(t) = +sin(t).
 *
 *
 * Это позволяет проверить второй порядок trapezoidal rule.
 */


/*
 * ============================================================
 * СОЗДАНИЕ МИНИМАЛЬНОЙ PEEC-МОДЕЛИ
 * ============================================================
 */
static int transient_test_build_model(
    Mesh *mesh,
    IncidenceMatrix *incidence,
    InductanceMatrix *L,
    PotentialMatrix *P,
    ResistanceMatrix *R,
    Matrix *PAT,
    int parallel_threads
)
{
    if (mesh == NULL ||
        incidence == NULL ||
        L == NULL ||
        P == NULL ||
        R == NULL ||
        PAT == NULL)
    {
        return -1;
    }

    mesh_init(mesh);
    incidence_init(incidence);
    inductance_init(L);
    potential_init(P);
    resistance_init(R);
    matrix_init(PAT);

    /*
     * --------------------------------------------------------
     * MESH
     * --------------------------------------------------------
     */
    mesh->n_nodes = 2;
    mesh->n_edges = 1;

    mesh->edges = malloc(
        2 * mesh->n_edges * sizeof(int)
    );

    if (mesh->edges == NULL)
    {
        return -1;
    }

    /*
     * Каноническая ориентация:
     *
     *     node 0 -> node 1.
     */
    mesh->edges[0] = 0;
    mesh->edges[1] = 1;


    /*
     * --------------------------------------------------------
     * INCIDENCE MATRIX
     * --------------------------------------------------------
     *
     *     A = [ -1  +1 ].
     */
    if (incidence_build(
            mesh,
            incidence
        ) != 0)
    {
        mesh_free(mesh);
        return -1;
    }


    /*
     * --------------------------------------------------------
     * L = [2]
     * --------------------------------------------------------
     */
    L->n = 1;

    L->data = calloc(
        1,
        sizeof(double)
    );

    if (L->data == NULL)
    {
        incidence_free(incidence);
        mesh_free(mesh);
        return -1;
    }

    L->data[0] = 2.0;


    /*
     * --------------------------------------------------------
     * P = I_2
     * --------------------------------------------------------
     */
    P->n = 2;

    P->data = calloc(
        4,
        sizeof(double)
    );

    if (P->data == NULL)
    {
        inductance_free(L);
        incidence_free(incidence);
        mesh_free(mesh);
        return -1;
    }

    P->data[0] = 1.0;
    P->data[3] = 1.0;


    /*
     * --------------------------------------------------------
     * R = 0
     * --------------------------------------------------------
     */
    R->n = 1;

    R->diagonal = calloc(
        1,
        sizeof(double)
    );

    if (R->diagonal == NULL)
    {
        potential_free(P);
        inductance_free(L);
        incidence_free(incidence);
        mesh_free(mesh);
        return -1;
    }

    R->diagonal[0] = 0.0;


    /*
     * --------------------------------------------------------
     * P A^T
     * --------------------------------------------------------
     */
    if (incidence_build_PAT(
            mesh,
            P,
            PAT,
            parallel_threads
        ) != 0)
    {
        matrix_free(PAT);
        resistance_free(R);
        potential_free(P);
        inductance_free(L);
        incidence_free(incidence);
        mesh_free(mesh);
        return -1;
    }

    /*
     * Ожидаем:
     *
     *     P A^T =
     *
     *     [ -1 ]
     *     [ +1 ].
     */
    if (fabs(PAT->data[0] + 1.0) > 1.0e-14 ||
        fabs(PAT->data[1] - 1.0) > 1.0e-14)
    {
        fprintf(
            stderr,
            "ERROR: transient test P A^T is incorrect.\n"
        );

        matrix_free(PAT);
        resistance_free(R);
        potential_free(P);
        inductance_free(L);
        incidence_free(incidence);
        mesh_free(mesh);

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ МИНИМАЛЬНОЙ МОДЕЛИ
 * ============================================================
 */
static void transient_test_free_model(
    Mesh *mesh,
    IncidenceMatrix *incidence,
    InductanceMatrix *L,
    PotentialMatrix *P,
    ResistanceMatrix *R,
    Matrix *PAT
)
{
    matrix_free(PAT);
    resistance_free(R);
    potential_free(P);
    inductance_free(L);
    incidence_free(incidence);
    mesh_free(mesh);
}


/*
 * ============================================================
 * STAGE-1 INTEGRATION HELPER
 * ============================================================
 */
static int transient_test_integrate_stage1(
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *incidence,
    const PotentialMatrix *P,
    const Matrix *PAT,
    double dt,
    double final_time,
    int parallel_threads,
    double *I_final,
    double *V0_final,
    double *V1_final
)
{
    if (L == NULL ||
        R == NULL ||
        incidence == NULL ||
        P == NULL ||
        PAT == NULL ||
        I_final == NULL ||
        V0_final == NULL ||
        V1_final == NULL ||
        !isfinite(dt) ||
        !isfinite(final_time) ||
        dt <= 0.0 ||
        final_time <= 0.0)
    {
        return -1;
    }

    TransientSystem system;

    transient_system_init(
        &system
    );

    if (transient_system_create_stage1(
            &system,
            L,
            R,
            incidence,
            P,
            PAT,
            dt
        ) != 0)
    {
        transient_system_free(&system);
        return -1;
    }

    if (transient_factorize(
            &system,
            parallel_threads
        ) != 0)
    {
        transient_system_free(&system);
        return -1;
    }

    /*
     * Начальное условие:
     *
     *     I = 1
     *     V = 0.
     */
    double I[1] = {
        1.0
    };

    double V[2] = {
        0.0,
        0.0
    };

    size_t n_steps =
        (size_t)llround(
            final_time / dt
        );

    if (n_steps == 0)
    {
        transient_system_free(&system);
        return -1;
    }

    /*
     * Проверяем, что конечное время представимо целым
     * числом выбранных шагов.
     */
    double integrated_time =
        (double)n_steps * dt;

    if (fabs(integrated_time - final_time) >
        1.0e-12 * fmax(1.0, fabs(final_time)))
    {
        transient_system_free(&system);
        return -1;
    }

    for (size_t step = 0; step < n_steps; ++step)
    {
        /*
         * Источников нет:
         *
         *     I_s = 0
         *     U_e = 0.
         */
        if (transient_build_rhs_stage1(
                &system,
                I,
                V,
                NULL,
                NULL,
                NULL,
                NULL
            ) != 0)
        {
            transient_system_free(&system);
            return -1;
        }

        if (transient_solve_step(
                &system
            ) != 0)
        {
            transient_system_free(&system);
            return -1;
        }

        if (transient_extract_stage1_solution(
                &system,
                I,
                V
            ) != 0)
        {
            transient_system_free(&system);
            return -1;
        }
    }

    *I_final = I[0];
    *V0_final = V[0];
    *V1_final = V[1];

    transient_system_free(
        &system
    );

    return 0;
}


/*
 * ============================================================
 * TEST 7A. STAGE 1 / ВТОРОЙ ПОРЯДОК
 * ============================================================
 */
static int test_transient_stage1(
    int parallel_threads
)
{
    Mesh mesh;
    IncidenceMatrix incidence;
    InductanceMatrix L;
    PotentialMatrix P;
    ResistanceMatrix R;
    Matrix PAT;

    if (transient_test_build_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT,
            parallel_threads
        ) != 0)
    {
        return -1;
    }

    printf("\n");
    printf("============================================================\n");
    printf("TRANSIENT STAGE-1 SELF-TEST\n");
    printf("============================================================\n");

    const double final_time = 1.0;

    double I_exact = cos(final_time);
    double V0_exact = -sin(final_time);
    double V1_exact = +sin(final_time);

    double I_dt = 0.0;
    double V0_dt = 0.0;
    double V1_dt = 0.0;

    double I_dt2 = 0.0;
    double V0_dt2 = 0.0;
    double V1_dt2 = 0.0;

    /*
     * Грубая сетка:
     *
     *     dt = 0.1.
     */
    if (transient_test_integrate_stage1(
            &L,
            &R,
            &incidence,
            &P,
            &PAT,
            0.1,
            final_time,
            parallel_threads,
            &I_dt,
            &V0_dt,
            &V1_dt
        ) != 0)
    {
        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    /*
     * В два раза меньший шаг:
     *
     *     dt = 0.05.
     */
    if (transient_test_integrate_stage1(
            &L,
            &R,
            &incidence,
            &P,
            &PAT,
            0.05,
            final_time,
            parallel_threads,
            &I_dt2,
            &V0_dt2,
            &V1_dt2
        ) != 0)
    {
        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    /*
     * Максимальная ошибка по всем неизвестным.
     */
    double error_dt =
        fmax(
            fabs(I_dt - I_exact),
            fmax(
                fabs(V0_dt - V0_exact),
                fabs(V1_dt - V1_exact)
            )
        );

    double error_dt2 =
        fmax(
            fabs(I_dt2 - I_exact),
            fmax(
                fabs(V0_dt2 - V0_exact),
                fabs(V1_dt2 - V1_exact)
            )
        );

    if (error_dt <= 0.0 ||
        error_dt2 <= 0.0)
    {
        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    /*
     * Для метода порядка p=2:
     *
     *         error(dt)
     *     ----------------
     *       error(dt/2)
     *
     * должно стремиться к:
     *
     *     2^p = 4.
     */
    double ratio =
        error_dt / error_dt2;

    printf(
        "Exact I(1)         : %.16e\n",
        I_exact
    );

    printf(
        "Computed dt=0.1    : %.16e\n",
        I_dt
    );

    printf(
        "Computed dt=0.05   : %.16e\n",
        I_dt2
    );

    printf("\n");

    printf(
        "Error dt           : %.9e\n",
        error_dt
    );

    printf(
        "Error dt/2         : %.9e\n",
        error_dt2
    );

    printf(
        "Error ratio        : %.9f\n",
        ratio
    );

    if (!isfinite(ratio) ||
        ratio < 3.5 ||
        ratio > 4.5)
    {
        fprintf(
            stderr,
            "ERROR: transient scheme does not show second-order convergence.\n"
        );

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    printf("Temporal order     : PASSED\n");
    printf("Stage-1 transient  : PASSED\n");
    printf("============================================================\n");

    transient_test_free_model(
        &mesh,
        &incidence,
        &L,
        &P,
        &R,
        &PAT
    );

    return 0;
}


/*
 * ============================================================
 * TEST 7B. STAGE 2 / SHUNT
 * ============================================================
 *
 * Используем ту же минимальную модель.
 *
 *
 * Подключаем внутренний шунт:
 *
 *     R_shunt = 2 Ohm
 *
 *     A = node 0
 *     B = node 1.
 *
 *
 * Положительное направление:
 *
 *     A -> B.
 *
 *
 * Тогда:
 *
 *     b_R =
 *
 *     [ -1 ]
 *     [ +1 ].
 *
 *
 * Начальное состояние:
 *
 *     I^0       = 1
 *     V_A^0     = 0
 *     V_B^0     = 0
 *     I_shunt^0 = 0.
 *
 *
 * При:
 *
 *     dt = 0.1
 *
 * после одного trapezoidal шага точное решение:
 *
 *     I^1       =  419 / 421
 *
 *     V_A^1     =  -40 / 421
 *
 *     V_B^1     =  +40 / 421
 *
 *     I_shunt^1 =  -40 / 421.
 */
static int test_transient_stage2(
    int parallel_threads
)
{
    Mesh mesh;
    IncidenceMatrix incidence;
    InductanceMatrix L;
    PotentialMatrix P;
    ResistanceMatrix R;
    Matrix PAT;

    if (transient_test_build_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT,
            parallel_threads
        ) != 0)
    {
        return -1;
    }

    printf("\n");
    printf("============================================================\n");
    printf("TRANSIENT STAGE-2 SELF-TEST\n");
    printf("============================================================\n");

    InternalCircuitShunt shunt;

    internal_circuit_shunt_default(
        &shunt
    );

    shunt.resistance = 2.0;
    shunt.node_a = 0;
    shunt.node_b = 1;
    shunt.nodes_mapped = 1;
    shunt.initialized = 1;

    TransientSystem system;

    transient_system_init(
        &system
    );

    if (transient_system_create_stage2(
            &system,
            &L,
            &R,
            &incidence,
            &P,
            &PAT,
            &shunt,
            0.1
        ) != 0)
    {
        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    /*
     * --------------------------------------------------------
     * b_R
     * --------------------------------------------------------
     */
    if (fabs(system.b_R[0] + 1.0) > 1.0e-14 ||
        fabs(system.b_R[1] - 1.0) > 1.0e-14)
    {
        fprintf(
            stderr,
            "ERROR: transient b_R sign is incorrect.\n"
        );

        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    /*
     * P = I:
     *
     *     P b_R = b_R.
     */
    if (fabs(system.P_bR[0] + 1.0) > 1.0e-14 ||
        fabs(system.P_bR[1] - 1.0) > 1.0e-14)
    {
        fprintf(
            stderr,
            "ERROR: transient P b_R is incorrect.\n"
        );

        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    if (transient_factorize(
            &system,
            parallel_threads
        ) != 0)
    {
        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    double I[1] = {
        1.0
    };

    double V[2] = {
        0.0,
        0.0
    };

    double I_shunt =
        0.0;

    /*
     * Внешних источников в unit-test нет.
     */
    if (transient_build_rhs_stage2(
            &system,
            I,
            V,
            I_shunt,
            NULL,
            NULL,
            NULL,
            NULL
        ) != 0)
    {
        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    if (transient_solve_step(
            &system
        ) != 0)
    {
        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    if (transient_extract_stage2_solution(
            &system,
            I,
            V,
            &I_shunt
        ) != 0)
    {
        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    const double I_expected =
        419.0 / 421.0;

    const double V0_expected =
        -40.0 / 421.0;

    const double V1_expected =
        +40.0 / 421.0;

    const double Ish_expected =
        -40.0 / 421.0;

    const double tolerance =
        1.0e-12;

    double max_error =
        fabs(
            I[0] - I_expected
        );

    double error =
        fabs(
            V[0] - V0_expected
        );

    if (error > max_error)
    {
        max_error = error;
    }

    error =
        fabs(
            V[1] - V1_expected
        );

    if (error > max_error)
    {
        max_error = error;
    }

    error =
        fabs(
            I_shunt - Ish_expected
        );

    if (error > max_error)
    {
        max_error = error;
    }

    double shunt_residual =
        internal_circuit_shunt_residual(
            &shunt,
            V,
            2,
            I_shunt
        );

    printf(
        "I                  : %.16e\n",
        I[0]
    );

    printf(
        "Expected I         : %.16e\n",
        I_expected
    );

    printf(
        "V_A                : %.16e V\n",
        V[0]
    );

    printf(
        "V_B                : %.16e V\n",
        V[1]
    );

    printf(
        "I_shunt            : %.16e A\n",
        I_shunt
    );

    printf(
        "Max solution error : %.9e\n",
        max_error
    );

    printf(
        "Shunt residual     : %.9e V\n",
        shunt_residual
    );

    if (max_error > tolerance ||
        !isfinite(shunt_residual) ||
        fabs(shunt_residual) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: transient Stage-2 self-test failed.\n"
        );

        transient_system_free(&system);

        transient_test_free_model(
            &mesh,
            &incidence,
            &L,
            &P,
            &R,
            &PAT
        );

        return -1;
    }

    printf("b_R sign           : PASSED\n");
    printf("P b_R              : PASSED\n");
    printf("Shunt equation     : PASSED\n");
    printf("Stage-2 transient  : PASSED\n");
    printf("============================================================\n");

    transient_system_free(
        &system
    );

    transient_test_free_model(
        &mesh,
        &incidence,
        &L,
        &P,
        &R,
        &PAT
    );

    return 0;
}


/*
 * ============================================================
 * TRANSIENT MASTER SELF-TEST
 * ============================================================
 */
static int test_transient(
    int parallel_threads
)
{
    if (test_transient_stage1(
            parallel_threads
        ) != 0)
    {
        return -1;
    }

    if (test_transient_stage2(
            parallel_threads
        ) != 0)
    {
        return -1;
    }

    printf("\n");
    printf("============================================================\n");
    printf("TRANSIENT PEEC SELF-TEST : PASSED\n");
    printf("============================================================\n");

    return 0;
}


/*
 * ============================================================
 * TEST 8. ПАДАЮЩЕЕ ПОЛЕ И ВЕТВЕВОЕ ВОЗБУЖДЕНИЕ
 * ============================================================
 */
static int test_scattering_excitation(
    int parallel_threads
)
{
    const double tolerance = 1.0e-12;

    Mesh mesh;
    IncidentField field;
    ScatteringOptions options;

    mesh_init(&mesh);
    incident_field_init(&field);

    scattering_options_default(
        &options,
        parallel_threads
    );

    mesh.n_edges = 2;

    mesh.edge_vectors = calloc(
        3 * mesh.n_edges,
        sizeof(double)
    );

    mesh.edge_centers = calloc(
        3 * mesh.n_edges,
        sizeof(double)
    );

    if (mesh.edge_vectors == NULL ||
        mesh.edge_centers == NULL)
    {
        mesh_free(&mesh);
        return -1;
    }

    /*
     * Edge 0:
     *
     *     l0 = (2,0,0).
     */
    mesh.edge_vectors[0] = 2.0;
    mesh.edge_vectors[1] = 0.0;
    mesh.edge_vectors[2] = 0.0;

    mesh.edge_centers[0] = 0.0;
    mesh.edge_centers[1] = 0.0;
    mesh.edge_centers[2] = 0.0;

    /*
     * Edge 1:
     *
     *     l1 = (0,3,0).
     */
    mesh.edge_vectors[3] = 0.0;
    mesh.edge_vectors[4] = 3.0;
    mesh.edge_vectors[5] = 0.0;

    mesh.edge_centers[3] = 0.0;
    mesh.edge_centers[4] = 0.0;
    mesh.edge_centers[5] = 0.0;

    printf("\n");
    printf("============================================================\n");
    printf("SCATTERING EXCITATION SELF-TEST\n");
    printf("============================================================\n");

    /*
     * Вертикальная поляризация:
     *
     * theta = 90 deg
     * phi   = 0 deg
     *
     * e_theta = (0,0,-1).
     *
     * Оба ребра лежат в xy -> U = 0.
     */
    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_COS,
            INCIDENT_POLARIZATION_VERTICAL,
            1.0,
            1.0e6,
            90.0,
            0.0,
            0.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    double U_time[2];

    if (scattering_build_time_excitation(
            &mesh,
            &field,
            &options,
            0.0,
            U_time
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_time[0]) > tolerance ||
        fabs(U_time[1]) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: vertical polarization scattering test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("Vertical polarization  : PASSED\n");

    /*
     * Горизонтальная:
     *
     * e_phi = (0,1,0)
     *
     * E_amp = 2 V/m
     *
     * edge1 = 3 m along y
     *
     * U1 = 6 V.
     */
    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_COS,
            INCIDENT_POLARIZATION_HORIZONTAL,
            2.0,
            1.0e6,
            90.0,
            0.0,
            0.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (scattering_build_time_excitation(
            &mesh,
            &field,
            &options,
            0.0,
            U_time
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_time[0]) > tolerance ||
        fabs(U_time[1] - 6.0) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: horizontal polarization scattering test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("Horizontal polarization: PASSED\n");

    /*
     * Проверка ориентации.
     */
    mesh.edge_vectors[4] = -3.0;

    if (scattering_build_time_excitation(
            &mesh,
            &field,
            &options,
            0.0,
            U_time
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_time[1] + 6.0) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: edge orientation scattering test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("Edge orientation sign  : PASSED\n");

    mesh.edge_vectors[4] = 3.0;

    /*
     * COS.
     */
    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_COS,
            INCIDENT_POLARIZATION_HORIZONTAL,
            2.0,
            1.0e6,
            90.0,
            0.0,
            0.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (scattering_build_time_excitation(
            &mesh,
            &field,
            &options,
            0.0,
            U_time
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_time[1] - 6.0) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: COS scattering test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("COS field               : PASSED\n");

    /*
     * SIN.
     */
    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_SIN,
            INCIDENT_POLARIZATION_HORIZONTAL,
            2.0,
            1.0e6,
            90.0,
            0.0,
            0.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (scattering_build_time_excitation(
            &mesh,
            &field,
            &options,
            0.0,
            U_time
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_time[1]) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: SIN scattering test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("SIN field               : PASSED\n");

    /*
     * Комплексное представление:
     *
     *     exp(j*w*t).
     */
    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_COMPLEX,
            INCIDENT_POLARIZATION_HORIZONTAL,
            2.0,
            1.0e6,
            90.0,
            0.0,
            0.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    Complex U_complex[2];

    if (scattering_build_harmonic_excitation(
            &mesh,
            &field,
            &options,
            U_complex
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_complex[1].re - 6.0) > tolerance ||
        fabs(U_complex[1].im) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: complex scattering excitation test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("Complex exp(j*w*t)      : PASSED\n");

    /*
     * Complex phase = pi/2.
     *
     * U = j*6 V.
     */
    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_COMPLEX,
            INCIDENT_POLARIZATION_HORIZONTAL,
            2.0,
            1.0e6,
            90.0,
            0.0,
            M_PI / 2.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (scattering_build_harmonic_excitation(
            &mesh,
            &field,
            &options,
            U_complex
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_complex[1].re) > tolerance ||
        fabs(U_complex[1].im - 6.0) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: complex scattering phase test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("Complex phase           : PASSED\n");

    /*
     * Импульс:
     *
     *     E(t) =
     *
     *     10 * (exp(-t) - exp(-2t)).
     */
    if (incident_field_configure_pulse(
            &field,
            INCIDENT_POLARIZATION_HORIZONTAL,
            10.0,
            1.0,
            2.0,
            90.0,
            0.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    double pulse_time = 1.0;

    if (scattering_build_time_excitation(
            &mesh,
            &field,
            &options,
            pulse_time,
            U_time
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    double expected_pulse =
        3.0 * 10.0 * (
            exp(-1.0)
            - exp(-2.0)
        );

    if (fabs(U_time[1] - expected_pulse) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: pulse scattering test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("Pulse field             : PASSED\n");

    /*
     * Причинность импульса.
     */
    mesh.edge_centers[3] = -299792458.0;
    mesh.edge_centers[4] = 0.0;
    mesh.edge_centers[5] = 0.0;

    if (scattering_build_time_excitation(
            &mesh,
            &field,
            &options,
            0.0,
            U_time
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }

    if (fabs(U_time[1]) > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: pulse causality scattering test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }

    printf("Pulse causality         : PASSED\n");
    /*
    * ========================================================
    * ТОЧНЫЙ ГАРМОНИЧЕСКИЙ ИНТЕГРАЛ ПО РЕБРУ
    * ========================================================
    *
    * Берем:
    *
    *     theta = 90 deg
    *     phi   = 45 deg.
    *
    *
    * Тогда:
    *
    *     k_hat =
    *
    *     (-1/sqrt(2), -1/sqrt(2), 0)
    *
    *
    * и:
    *
    *     e_phi =
    *
    *     (-1/sqrt(2), +1/sqrt(2), 0).
    *
    *
    * Ребро:
    *
    *     l = (0,3,0).
    *
    *
    * Значит:
    *
    *     E0.l =
    *
    *     2 * 3 / sqrt(2).
    *
    *
    * Для центра в начале координат:
    *
    *     E_hat(r_c) = E0.
    *
    *
    * Фазовый параметр:
    *
    *     q = k0 * k_hat.l.
    *
    *
    * Точное значение:
    *
    *     U =
    *
    *     (E0.l) * sinc(q/2).
    */
    mesh.edge_centers[3] = 0.0;
    mesh.edge_centers[4] = 0.0;
    mesh.edge_centers[5] = 0.0;

    mesh.edge_vectors[3] = 0.0;
    mesh.edge_vectors[4] = 3.0;
    mesh.edge_vectors[5] = 0.0;


    double test_frequency =
        50.0e6;

    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_COMPLEX,
            INCIDENT_POLARIZATION_HORIZONTAL,
            2.0,
            test_frequency,
            90.0,
            45.0,
            0.0
        ) != 0)
    {
        mesh_free(&mesh);
        return -1;
    }


    options.integration =
        SCATTERING_EDGE_ANALYTIC_PLANE_WAVE;


    if (scattering_build_harmonic_excitation(
            &mesh,
            &field,
            &options,
            U_complex
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: analytic edge integral test failed.\n"
        );

        mesh_free(&mesh);
        return -1;
    }


    /*
    * Строим ожидаемое значение независимо
    * от scattering.c.
    */
    double inv_sqrt2 =
        1.0 / sqrt(2.0);

    double e_dot_l =
        2.0 * 3.0 * inv_sqrt2;

    double k_dot_l =
        -3.0 * inv_sqrt2;

    double q =
        field.wave_number
        * k_dot_l;

    double x =
        0.5 * q;

    double sinc_value;

    if (fabs(x) < 1.0e-12)
    {
        sinc_value = 1.0;
    }
    else
    {
        sinc_value = sin(x) / x;
    }

    double expected_exact =
        e_dot_l * sinc_value;


    /*
    * Центр находится в начале координат, phase=0,
    * поэтому ожидаем чисто действительный фазор.
    */
    if (fabs(U_complex[1].re - expected_exact) > 1.0e-11 ||
        fabs(U_complex[1].im) > 1.0e-11)
    {
        fprintf(
            stderr,
            "ERROR: analytic plane-wave edge integral is incorrect.\n"
        );

        printf(
            "Computed U1 : %.16e %+.16ej V\n",
            U_complex[1].re,
            U_complex[1].im
        );

        printf(
            "Expected U1 : %.16e +0j V\n",
            expected_exact
        );

        mesh_free(&mesh);
        return -1;
    }

    printf(
        "Analytic edge integral : PASSED\n"
    );


    /*
    * Возвращаем MIDPOINT, чтобы последующие временные
    * тесты не зависели от гармонического режима.
    */
    options.integration =
        SCATTERING_EDGE_MIDPOINT;
    printf("============================================================\n");
    printf("SCATTERING EXCITATION TEST : PASSED\n");
    printf("============================================================\n");

    mesh_free(&mesh);

    return 0;
}

#endif


/*
 * ============================================================
 * MAIN
 * ============================================================
 */
int main(int argc, char **argv)
{
    /*
     * ========================================================
     * 1. КОНФИГУРАЦИЯ
     * ========================================================
     */
    PeecConfig config;

    config_set_defaults(&config);

    int parse_result =
        config_parse_cli(
            &config,
            argc,
            argv
        );

    if (parse_result > 0)
    {
        return EXIT_SUCCESS;
    }

    if (parse_result < 0)
    {
        config_print_help(argv[0]);
        return EXIT_FAILURE;
    }


    /*
     * ========================================================
     * 2. OPENMP
     * ========================================================
     */
    omp_set_dynamic(0);
    omp_set_num_threads(config.parallel_threads);


#ifdef TEST

    /*
     * ========================================================
     * 3. NUMERICAL SELF-TESTS
     * ========================================================
     */
    if (test_lu(config.parallel_threads) != 0)
    {
        fprintf(
            stderr,
            "ERROR: real LU self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    if (test_complex_arithmetic() != 0)
    {
        fprintf(
            stderr,
            "ERROR: complex arithmetic self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    if (test_complex_matrix(config.parallel_threads) != 0)
    {
        fprintf(
            stderr,
            "ERROR: complex matrix self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    if (test_complex_lu(config.parallel_threads) != 0)
    {
        fprintf(
            stderr,
            "ERROR: complex LU self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    if (test_resistance(config.parallel_threads) != 0)
    {
        fprintf(
            stderr,
            "ERROR: resistance self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    if (test_incidence_PAT(config.parallel_threads) != 0)
    {
        fprintf(
            stderr,
            "ERROR: P A^T self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    if (test_transient(
            config.parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: transient PEEC self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    if (test_scattering_excitation(
            config.parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: scattering excitation self-test failed.\n"
        );

        return EXIT_FAILURE;
    }

    printf("\n");
    printf("============================================================\n");
    printf("ALL NUMERICAL SELF-TESTS PASSED\n");
    printf("============================================================\n");

#endif


    /*
     * ========================================================
     * 4. ИНФОРМАЦИЯ О ЗАПУСКЕ
     * ========================================================
     */
    printf("\n");
    printf("============================================================\n");
    printf("PEEC / PARTIAL EQUIVALENT CIRCUIT\n");
    printf("============================================================\n");

    printf(
        "Task        : %s\n",
        task_name(config.task)
    );

    printf(
        "Physics     : %s\n",
        physics_name(config.physics)
    );

    printf(
        "Mesh        : %s\n",
        config.mesh_file
    );

    printf(
        "Parallel    : %d OpenMP thread(s)\n",
        config.parallel_threads
    );

    printf(
        "OMP max     : %d\n",
        omp_get_max_threads()
    );

#ifdef TEST
    printf("Build       : TEST\n");
#endif

    printf("============================================================\n");


    /*
     * ========================================================
     * 5. СЕТКА
     * ========================================================
     */
    Mesh mesh;

    mesh_init(&mesh);

    printf("\nReading Gmsh mesh ...\n");

    if (mesh_load_gmsh(
            config.mesh_file,
            &mesh
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: mesh loading failed.\n"
        );

        mesh_free(&mesh);

        return EXIT_FAILURE;
    }


    /*
     * ========================================================
     * 6. ТОПОЛОГИЯ
     * ========================================================
     */
    printf("Building topology ...\n");

    if (mesh_build_topology(
            &mesh
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: topology construction failed.\n"
        );

        mesh_free(&mesh);

        return EXIT_FAILURE;
    }

    mesh_print_info(&mesh);


#ifdef TEST

    long long euler =
        (long long)mesh.n_nodes
        - (long long)mesh.n_edges
        + (long long)mesh.n_quads;

    printf(
        "Euler characteristic : %lld\n",
        euler
    );

#endif


    /*
     * ========================================================
     * 7. ДУАЛЬНАЯ СЕТКА
     * ========================================================
     */
    DualMesh dual;

    dual_mesh_init(&dual);

    printf("\nBuilding dual mesh ...\n");

    if (dual_mesh_build(
            &mesh,
            &dual
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: dual mesh construction failed.\n"
        );

        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return EXIT_FAILURE;
    }

    dual_mesh_print_info(&dual);


    /*
     * ========================================================
     * 8. МАТРИЦА ИНЦИДЕНТНОСТИ A
     * ========================================================
     *
     * Размер:
     *
     *     A : Ne x Nv.
     */
    IncidenceMatrix incidence;

    incidence_init(&incidence);

    printf("\nBuilding incidence matrix ...\n");

    if (incidence_build(
            &mesh,
            &incidence
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: incidence matrix construction failed.\n"
        );

        incidence_free(&incidence);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return EXIT_FAILURE;
    }


#ifdef TEST

    incidence_print_info(&incidence);

    double incidence_error =
        incidence_check_row_sum(&incidence);

    printf(
        "Max |A * 1| : %.3e\n",
        incidence_error
    );

    if (incidence_error > 1.0e-14)
    {
        fprintf(
            stderr,
            "ERROR: incidence matrix consistency check failed.\n"
        );

        incidence_free(&incidence);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return EXIT_FAILURE;
    }

#endif


    /*
     * ========================================================
     * 9. МАТРИЦА L
     * ========================================================
     */
    InductanceOptions inductance_options;

    inductance_options_default(
        &inductance_options,
        config.parallel_threads
    );

    InductanceMatrix L;

    inductance_init(&L);

    printf("\nComputing inductance matrix L ...\n");

    double inductance_start =
        omp_get_wtime();

    if (inductance_compute_matrix(
            &mesh,
            &dual,
            &inductance_options,
            &L
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: inductance matrix calculation failed.\n"
        );

        inductance_free(&L);
        incidence_free(&incidence);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return EXIT_FAILURE;
    }

    double inductance_elapsed =
        omp_get_wtime()
        - inductance_start;

#ifdef TEST
    inductance_print_info(&L);
#endif

    printf(
        "L computation time : %.6f s\n",
        inductance_elapsed
    );


    /*
     * ========================================================
     * 10. МАТРИЦА P
     * ========================================================
     */
    PotentialOptions potential_options;

    potential_options_default(
        &potential_options,
        config.parallel_threads
    );

    PotentialMatrix P;

    potential_init(&P);

    printf("\nComputing potential matrix P ...\n");

    double potential_start =
        omp_get_wtime();

    if (potential_compute_matrix(
            &mesh,
            &dual,
            &potential_options,
            &P
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: potential matrix calculation failed.\n"
        );

        potential_free(&P);
        inductance_free(&L);
        incidence_free(&incidence);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return EXIT_FAILURE;
    }

    double potential_elapsed =
        omp_get_wtime()
        - potential_start;

#ifdef TEST
    potential_print_info(&P);
#endif

    printf(
        "P computation time : %.6f s\n",
        potential_elapsed
    );


    /*
     * ========================================================
     * 11. P A^T
     * ========================================================
     *
     * Размер:
     *
     *     Nv x Ne.
     *
     *
     * Используется непосредственно в классической
     * PEEC-системе:
     *
     *     dV/dt - P A^T I = P I_s.
     */
    Matrix PAT;

    matrix_init(&PAT);

    printf("\nBuilding P A^T ...\n");

    double pat_start =
        omp_get_wtime();

    if (incidence_build_PAT(
            &mesh,
            &P,
            &PAT,
            config.parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: failed to build P A^T.\n"
        );

        matrix_free(&PAT);
        potential_free(&P);
        inductance_free(&L);
        incidence_free(&incidence);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return EXIT_FAILURE;
    }

    double pat_elapsed =
        omp_get_wtime()
        - pat_start;

    printf(
        "P A^T computation time : %.6f s\n",
        pat_elapsed
    );


    /*
     * ========================================================
     * 12. СОПРОТИВЛЕНИЕ R
     * ========================================================
     */
    ResistanceOptions resistance_options;

    resistance_options_default(
        &resistance_options,
        config.parallel_threads
    );

    /*
     * Пока основная задача решается для PEC:
     *
     *     R = 0.
     */
    resistance_options.model =
        PEEC_CONDUCTOR_PEC;

    ResistanceMatrix R;

    resistance_init(&R);

    printf("\nBuilding resistance matrix R ...\n");

    double resistance_start =
        omp_get_wtime();

    if (resistance_compute_matrix(
            &mesh,
            &dual,
            &resistance_options,
            &R
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: resistance matrix calculation failed.\n"
        );

        resistance_free(&R);
        matrix_free(&PAT);
        potential_free(&P);
        inductance_free(&L);
        incidence_free(&incidence);
        dual_mesh_free(&dual);
        mesh_free(&mesh);

        return EXIT_FAILURE;
    }

    double resistance_elapsed =
        omp_get_wtime()
        - resistance_start;

#ifdef TEST
    resistance_print_info(
        &R,
        &resistance_options
    );
#endif

    printf(
        "R computation time : %.6f s\n",
        resistance_elapsed
    );


    /*
     * ========================================================
     * 13. СВОДКА ПО PEEC-МАТРИЦАМ
     * ========================================================
     */
    printf("\n");
    printf("PEEC matrix assembly summary\n");
    printf("----------------------------------------\n");

    printf(
        "L time       : %.6f s\n",
        inductance_elapsed
    );

    printf(
        "P time       : %.6f s\n",
        potential_elapsed
    );

    printf(
        "P A^T time   : %.6f s\n",
        pat_elapsed
    );

    printf(
        "R time       : %.6f s\n",
        resistance_elapsed
    );

    printf(
        "Total        : %.6f s\n",
        inductance_elapsed
        + potential_elapsed
        + pat_elapsed
        + resistance_elapsed
    );


    printf("\n");
    printf("Classical PEEC system\n");
    printf("----------------------------------------\n");

    printf(
        "Edge unknowns I : %zu\n",
        mesh.n_edges
    );

    printf(
        "Node unknowns V : %zu\n",
        mesh.n_nodes
    );

    printf(
        "Full system size: %zu x %zu\n",
        mesh.n_edges + mesh.n_nodes,
        mesh.n_edges + mesh.n_nodes
    );


    /*
     * ========================================================
     * 14. ВЫБОР ФИЗИЧЕСКОЙ ЗАДАЧИ
     * ========================================================
     */
    switch (config.task)
    {
        /*
         * ----------------------------------------------------
         * MESH INFO
         * ----------------------------------------------------
         */
        case PEEC_TASK_MESH_INFO:
        {
            break;
        }


        /*
         * ----------------------------------------------------
         * ГАРМОНИЧЕСКОЕ РАССЕЯНИЕ
         * ----------------------------------------------------
         *
         * Используем полную классическую систему:
         *
         *     [ R+j*w*L      A     ] [ I ]   [ -U_inc ]
         *     [                    ] [   ] = [        ]
         *     [ -P*A^T      j*w*Iv ] [ V ]   [   0    ].
         *
         *
         * Пока параметры падающей волны заданы здесь.
         * Позже вынесем их в config/case.
         */
        case PEEC_TASK_SCATTERING:
        {
            printf("\n");
            printf("============================================================\n");
            printf("HARMONIC SCATTERING\n");
            printf("============================================================\n");


            /*
             * =================================================
             * ПАРАМЕТРЫ ПАДАЮЩЕЙ ВОЛНЫ
             * =================================================
             *
             * Сейчас:
             *
             *     E_amp = 1 V/m
             *     f     = 300 MHz
             *     theta = 90 deg
             *     phi   = 0 deg
             *
             * Горизонтальная поляризация:
             *
             *     E || e_phi.
             *
             * Комплексное представление:
             *
             *     exp(j*w*t).
             */
            IncidentField field;

            incident_field_init(&field);

            double frequency =
                300.0e6;

            if (incident_field_configure_harmonic(
                    &field,
                    INCIDENT_WAVE_COMPLEX,
                    INCIDENT_POLARIZATION_HORIZONTAL,
                    1.0,
                    frequency,
                    90.0,
                    45.0,
                    0.0
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: cannot configure incident field.\n"
                );

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }

            incident_field_print_info(
                &field
            );


            /*
             * =================================================
             * U_inc ПО РЕБРАМ
             * =================================================
             */
            ScatteringOptions scattering_options;

            scattering_options_default(
                &scattering_options,
                config.parallel_threads
            );
            /*
            * Для монохроматической плоской волны используем
            * точный аналитический интеграл вдоль ребра:
            *
            *     U_e = integral E_inc . dl.
            */
            scattering_options.integration =SCATTERING_EDGE_ANALYTIC_PLANE_WAVE;
            Complex *U_inc =
                malloc(
                    mesh.n_edges
                    * sizeof(Complex)
                );

            if (U_inc == NULL)
            {
                fprintf(
                    stderr,
                    "ERROR: cannot allocate harmonic incident excitation.\n"
                );

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }

            if (scattering_build_harmonic_excitation(
                    &mesh,
                    &field,
                    &scattering_options,
                    U_inc
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: cannot build harmonic scattering excitation.\n"
                );

                free(U_inc);

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }

            scattering_print_harmonic_excitation_info(
                U_inc,
                mesh.n_edges
            );


            /*
             * =================================================
             * ПОЛНАЯ ГАРМОНИЧЕСКАЯ PEEC-СИСТЕМА
             * =================================================
             */
            HarmonicSystem harmonic_system;

            harmonic_system_init(
                &harmonic_system
            );

            printf(
                "\nBuilding full harmonic PEEC system ...\n"
            );

            double harmonic_build_start =
                omp_get_wtime();

            if (harmonic_build_system(
                    &L,
                    &R,
                    &incidence,
                    &PAT,
                    frequency,
                    config.parallel_threads,
                    &harmonic_system
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: harmonic PEEC system construction failed.\n"
                );

                harmonic_system_free(
                    &harmonic_system
                );

                free(U_inc);

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }

            double harmonic_build_elapsed =
                omp_get_wtime()
                - harmonic_build_start;

            harmonic_print_info(
                &harmonic_system
            );

            printf(
                "Harmonic matrix build time : %.6f s\n",
                harmonic_build_elapsed
            );


            /*
             * =================================================
             * LU-ФАКТОРИЗАЦИЯ
             * =================================================
             */
            printf(
                "\nFactorizing full harmonic PEEC system ...\n"
            );

            double harmonic_lu_start =
                omp_get_wtime();

            if (harmonic_factorize(
                    &harmonic_system,
                    config.parallel_threads
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: harmonic PEEC factorization failed.\n"
                );

                harmonic_system_free(
                    &harmonic_system
                );

                free(U_inc);

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }

            double harmonic_lu_elapsed =
                omp_get_wtime()
                - harmonic_lu_start;

            printf(
                "Harmonic LU time           : %.6f s\n",
                harmonic_lu_elapsed
            );


            /*
             * =================================================
             * ВОЗБУЖДЕНИЕ
             * =================================================
             *
             * Для задачи рассеяния:
             *
             *     U_e = U_inc,
             *
             *     I_s = 0.
             */
            HarmonicExcitation excitation;

            excitation.n_edges =
                mesh.n_edges;

            excitation.n_nodes =
                mesh.n_nodes;

            excitation.edge_voltage =
                U_inc;

            excitation.node_current =
                NULL;


            /*
             * =================================================
             * РЕШЕНИЕ
             * =================================================
             */
            HarmonicSolution solution;

            harmonic_solution_init(
                &solution
            );

            printf(
                "\nSolving harmonic PEEC system ...\n"
            );

            double harmonic_solve_start =
                omp_get_wtime();

            if (harmonic_solve(
                    &harmonic_system,
                    &P,
                    &excitation,
                    &solution,
                    config.parallel_threads
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: harmonic PEEC solve failed.\n"
                );

                harmonic_solution_free(
                    &solution
                );

                harmonic_system_free(
                    &harmonic_system
                );

                free(U_inc);

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }

            double harmonic_solve_elapsed =
                omp_get_wtime()
                - harmonic_solve_start;


            /*
             * =================================================
             * НЕВЯЗКА ПОЛНОЙ PEEC-СИСТЕМЫ
             * =================================================
             */
            HarmonicResidual residual;

            if (harmonic_compute_residual(
                    &harmonic_system,
                    &P,
                    &excitation,
                    &solution,
                    config.parallel_threads,
                    &residual
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: harmonic PEEC residual calculation failed.\n"
                );

                harmonic_solution_free(
                    &solution
                );

                harmonic_system_free(
                    &harmonic_system
                );

                free(U_inc);

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }


            /*
             * =================================================
             * ДИАГНОСТИКА РЕШЕНИЯ
             * =================================================
             */
            double max_current =
                complex_vector_max_abs(
                    solution.edge_current,
                    solution.n_edges
                );

            double max_voltage =
                complex_vector_max_abs(
                    solution.node_voltage,
                    solution.n_nodes
                );


            printf("\n");
            printf("========================================\n");
            printf("HARMONIC PEEC SOLUTION\n");
            printf("========================================\n");

            printf(
                "Frequency          : %.9e Hz\n",
                frequency
            );

            printf(
                "Maximum |I_e|      : %.9e A\n",
                max_current
            );

            printf(
                "Maximum |V_j|      : %.9e V\n",
                max_voltage
            );

            printf("\n");

            /*
             * Эта величина объединяет два блока с разными
             * физическими размерностями, поэтому используется
             * только как вспомогательная диагностика.
             */
            printf(
                "Absolute residual  : %.9e\n",
                residual.absolute_residual
            );

            /*
             * Также только дополнительная диагностика.
             */
            printf(
                "RHS-relative       : %.9e\n",
                residual.rhs_relative_residual
            );

            /*
             * Главный численный критерий качества решения.
             */
            printf(
                "Backward error     : %.9e\n",
                residual.backward_error
            );

            printf("\n");

            /*
             * Верхний PEEC-блок имеет размерность В.
             */
            printf(
                "Edge residual      : %.9e V\n",
                residual.edge_absolute_residual
            );

            /*
             * Нижний блок:
             *
             *     -P A^T I + j*w V
             *
             * имеет размерность В/с.
             */
            printf(
                "Node residual      : %.9e V/s\n",
                residual.node_absolute_residual
            );

            printf(
                "Edge back. error   : %.9e\n",
                residual.edge_backward_error
            );

            printf(
                "Node back. error   : %.9e\n",
                residual.node_backward_error
            );

            printf("\n");

            printf(
                "Solve time         : %.6f s\n",
                harmonic_solve_elapsed
            );

            printf(
                "Matrix build       : %.6f s\n",
                harmonic_build_elapsed
            );

            printf(
                "LU factorization   : %.6f s\n",
                harmonic_lu_elapsed
            );

            printf("========================================\n");


#ifdef TEST

            /*
             * =================================================
             * ПРОВЕРКА КАЧЕСТВА РЕШЕНИЯ
             * =================================================
             *
             * Не используем:
             *
             *     ||Zx-b|| / ||b||
             *
             * как основной критерий, потому что масштабы
             * двух PEEC-блоков сильно различаются.
             *
             *
             * Используем компонентно нормированную:
             *
             *                    |r_i|
             *     eta = max -----------------------------
             *            i  sum_j |Z_ij||x_j| + |b_i|.
             */
            if (!isfinite(residual.backward_error) ||
                residual.backward_error > 1.0e-10)
            {
                fprintf(
                    stderr,
                    "ERROR: harmonic PEEC backward error is too large: %.9e\n",
                    residual.backward_error
                );

                harmonic_solution_free(
                    &solution
                );

                harmonic_system_free(
                    &harmonic_system
                );

                free(U_inc);

                resistance_free(&R);
                matrix_free(&PAT);
                potential_free(&P);
                inductance_free(&L);
                incidence_free(&incidence);
                dual_mesh_free(&dual);
                mesh_free(&mesh);

                return EXIT_FAILURE;
            }

            printf(
                "Harmonic residual test : PASSED\n"
            );

#endif


            /*
             * =================================================
             * ОСВОБОЖДЕНИЕ ЛОКАЛЬНЫХ ДАННЫХ
             * =================================================
             */
            harmonic_solution_free(
                &solution
            );

            harmonic_system_free(
                &harmonic_system
            );

            free(U_inc);

            break;
        }


        /*
         * ----------------------------------------------------
         * МОЛНИЯ
         * ----------------------------------------------------
         *
         * Следующий solver:
         *
         *     transient.c.
         *
         * Для прямого удара молнии возбуждение будет
         * задаваться через узловой ток:
         *
         *     I_s(t).
         */
        case PEEC_TASK_LIGHTNING:
        {
            printf("\n");
            printf(
                "Lightning transient PEEC solver is not implemented yet.\n"
            );

            break;
        }


        case PEEC_TASK_NONE:
        default:
        {
            fprintf(
                stderr,
                "ERROR: invalid task.\n"
            );

            resistance_free(&R);
            matrix_free(&PAT);
            potential_free(&P);
            inductance_free(&L);
            incidence_free(&incidence);
            dual_mesh_free(&dual);
            mesh_free(&mesh);

            return EXIT_FAILURE;
        }
    }


    /*
     * ========================================================
     * 15. ОСВОБОЖДЕНИЕ PEEC-МАТРИЦ
     * ========================================================
     */
    resistance_free(&R);
    matrix_free(&PAT);
    potential_free(&P);
    inductance_free(&L);
    incidence_free(&incidence);
    dual_mesh_free(&dual);
    mesh_free(&mesh);


    printf("\n");
    printf("Done.\n");

    return EXIT_SUCCESS;
}