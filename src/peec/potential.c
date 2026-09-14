#include "potential.h"
#include "quadrature.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * ============================================================
 * ФИЗИЧЕСКАЯ КОНСТАНТА
 * ============================================================
 *
 * Электрическая постоянная:
 *
 *     eps0 = 8.8541878128e-12 Ф/м.
 *
 * Коэффициент:
 *
 *          1
 *     -----------
 *     4*pi*eps0
 *
 * имеет размерность:
 *
 *     В*м/Кл.
 */
static const double INV_4PI_EPS0 = 8.987551792261171e9;


/*
 * ============================================================
 * ПАРАМЕТРЫ ПО УМОЛЧАНИЮ
 * ============================================================
 */
void potential_options_default(
    PotentialOptions *options,
    int parallel_threads
)
{
    if (options == NULL)
    {
        return;
    }

    /*
     * Используем те же уровни точности, что и для L.
     */
    options->far_order = 4;
    options->near_order = 10;
    options->touching_order = 12;
    options->self_order = 8;

    options->near_factor = 0.5;
    options->parallel_threads = parallel_threads;
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void potential_init(PotentialMatrix *P)
{
    if (P == NULL)
    {
        return;
    }

    memset(P, 0, sizeof(*P));
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ ПАМЯТИ
 * ============================================================
 */
void potential_free(PotentialMatrix *P)
{
    if (P == NULL)
    {
        return;
    }

    free(P->data);

    potential_init(P);
}


/*
 * ============================================================
 * ИНТЕГРАЛ ДЛЯ ДВУХ QUADPATCH
 * ============================================================
 */
static int compute_patch_pair_kernel(
    const QuadPatch *patch_a,
    const QuadPatch *patch_b,
    const PotentialOptions *options,
    double *value
)
{
    if (patch_a == NULL ||
        patch_b == NULL ||
        options == NULL ||
        value == NULL)
    {
        return -1;
    }

    *value = quad_patch_kernel_1_over_r(
        patch_a,
        patch_b,
        options->far_order,
        options->near_order,
        options->touching_order,
        options->self_order,
        options->near_factor
    );

    if (!isfinite(*value))
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ОДИН КОЭФФИЦИЕНТ P_ij
 * ============================================================
 */
int potential_compute_pair(
    const Mesh *mesh,
    const DualMesh *dual,
    size_t node_i,
    size_t node_j,
    const PotentialOptions *options,
    double *value
)
{
    if (mesh == NULL ||
        dual == NULL ||
        options == NULL ||
        value == NULL)
    {
        return -1;
    }

    if (node_i >= mesh->n_nodes || node_j >= mesh->n_nodes)
    {
        return -1;
    }

    /*
     * Площади узловых дуальных областей:
     *
     *     S_i^v
     *     S_j^v.
     */
    double area_i = dual->node_region_areas[node_i];
    double area_j = dual->node_region_areas[node_j];

    if (area_i <= 0.0 || area_j <= 0.0)
    {
        return -1;
    }

    /*
     * Узловая область может состоять из нескольких
     * quarter-patch.
     */
    const PatchList *region_i =
    &dual->node_patches[node_i];

    const PatchList *region_j =
    &dual->node_patches[node_j];

    if (region_i->count == 0 || region_j->count == 0)
    {
        return -1;
    }

    /*
     * Геометрический интеграл:
     *
     *                  /        /
     *                 |        |       1
     *     G_ij =      |        |    ------- dS dS'.
     *                 |        |    |r-r'|
     *                /Pi_i    /Pi_j
     */
    double integral = 0.0;

    for (size_t a = 0; a < region_i->count; ++a)
    {
        const QuadPatch *patch_a =
        &region_i->patches[a];

        for (size_t b = 0; b < region_j->count; ++b)
        {
            const QuadPatch *patch_b =
            &region_j->patches[b];

            double patch_value = 0.0;

            if (compute_patch_pair_kernel(
                patch_a,
                patch_b,
                options,
                &patch_value
            ) != 0)
            {
                return -1;
            }

            integral += patch_value;
        }
    }

    /*
     * Полная формула:
     *
     *              1
     *     P_ij = --------- *
     *            4*pi*eps0
     *
     *              G_ij
     *          -----------
     *          S_i^v S_j^v.
     *
     * Размерность:
     *
     *     В/Кл = 1/Ф.
     */
    *value =
    INV_4PI_EPS0
    * integral
    / (area_i * area_j);

    if (!isfinite(*value))
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ПОЛНАЯ МАТРИЦА P
 * ============================================================
 */
int potential_compute_matrix(
    const Mesh *mesh,
    const DualMesh *dual,
    const PotentialOptions *options,
    PotentialMatrix *P
)
{
    if (mesh == NULL ||
        dual == NULL ||
        options == NULL ||
        P == NULL)
    {
        return -1;
    }

    if (mesh->n_nodes == 0)
    {
        return -1;
    }

    potential_init(P);

    size_t nv = mesh->n_nodes;

    P->n = nv;

    /*
     * Плотная матрица:
     *
     *     Nv x Nv.
     */
    P->data = calloc(
        nv * nv,
        sizeof(double)
    );

    if (P->data == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate potential matrix.\n"
        );

        return -1;
    }

    int error_code = 0;


    /*
     * ========================================================
     * OPENMP
     * ========================================================
     *
     * Используем симметрию:
     *
     *     P_ij = P_ji.
     *
     * schedule(dynamic) полезен, потому что различные пары
     * требуют разного порядка квадратуры.
     */
    #pragma omp parallel for schedule(dynamic) num_threads(options->parallel_threads)
    for (long long i_signed = 0; i_signed < (long long)nv; ++i_signed)
    {
        size_t i = (size_t)i_signed;

        int current_error = 0;

        #pragma omp atomic read
        current_error = error_code;

        if (current_error != 0)
        {
            continue;
        }

        for (size_t j = i; j < nv; ++j)
        {
            #pragma omp atomic read
            current_error = error_code;

            if (current_error != 0)
            {
                break;
            }

            double value = 0.0;

            int status = potential_compute_pair(
                mesh,
                dual,
                i,
                j,
                options,
                &value
            );

            if (status != 0)
            {
                #pragma omp critical(potential_error)
                {
                    if (error_code == 0)
                    {
                        error_code = status;

                        fprintf(
                            stderr,
                            "\nERROR: potential calculation failed "
                            "for nodes %zu, %zu.\n",
                            i,
                            j
                        );
                    }
                }

                break;
            }

            /*
             * Симметричное заполнение.
             */
            P->data[i * nv + j] = value;
            P->data[j * nv + i] = value;
        }
    }

    if (error_code != 0)
    {
        potential_free(P);

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ДИАГНОСТИКА МАТРИЦЫ P
 * ============================================================
 */
void potential_print_info(const PotentialMatrix *P)
{
    if (P == NULL || P->data == NULL || P->n == 0)
    {
        return;
    }

    double min_value = P->data[0];
    double max_value = P->data[0];

    double min_diagonal = HUGE_VAL;
    double max_diagonal = -HUGE_VAL;

    double max_symmetry_error = 0.0;

    for (size_t i = 0; i < P->n; ++i)
    {
        double diagonal = P->data[i * P->n + i];

        if (diagonal < min_diagonal)
        {
            min_diagonal = diagonal;
        }

        if (diagonal > max_diagonal)
        {
            max_diagonal = diagonal;
        }

        for (size_t j = 0; j < P->n; ++j)
        {
            double value = P->data[i * P->n + j];

            if (value < min_value)
            {
                min_value = value;
            }

            if (value > max_value)
            {
                max_value = value;
            }

            double symmetry_error = fabs(
                value - P->data[j * P->n + i]
            );

            if (symmetry_error > max_symmetry_error)
            {
                max_symmetry_error = symmetry_error;
            }
        }
    }

    double memory_mib =
    ((double)P->n
    * (double)P->n
    * (double)sizeof(double))
    / (1024.0 * 1024.0);

    printf("\n");
    printf("========================================\n");
    printf("POTENTIAL MATRIX\n");
    printf("========================================\n");
    printf("Size               : %zu x %zu\n", P->n, P->n);
    printf("Memory             : %.3f MiB\n", memory_mib);
    printf("Minimum P          : %.9e 1/F\n", min_value);
    printf("Maximum P          : %.9e 1/F\n", max_value);
    printf("Minimum diagonal   : %.9e 1/F\n", min_diagonal);
    printf("Maximum diagonal   : %.9e 1/F\n", max_diagonal);
    printf("Max symmetry error : %.3e 1/F\n", max_symmetry_error);
    printf("========================================\n");
}


/*
 * ============================================================
 * F_j = 1 / P_jj
 * ============================================================
 */
int potential_build_self_capacitance(
    const PotentialMatrix *P,
    double *F
)
{
    if (P == NULL ||
        P->data == NULL ||
        P->n == 0 ||
        F == NULL)
    {
        return -1;
    }

    for (size_t j = 0; j < P->n; ++j)
    {
        double diagonal =
        P->data[j * P->n + j];

        if (!isfinite(diagonal) ||
            diagonal <= 0.0)
        {
            fprintf(
                stderr,
                "ERROR: invalid diagonal P[%zu,%zu] in "
                "potential_build_self_capacitance().\n",
                    j,
                    j
            );

            return -1;
        }

        F[j] =
        1.0 / diagonal;
    }

    return 0;
}


/*
 * ============================================================
 * q = F V_c
 * ============================================================
 */
int potential_self_voltage_to_charge(
    const PotentialMatrix *P,
    const double *V_c,
    double *charge
)
{
    if (P == NULL ||
        P->data == NULL ||
        P->n == 0 ||
        V_c == NULL ||
        charge == NULL)
    {
        return -1;
    }

    for (size_t j = 0; j < P->n; ++j)
    {
        double diagonal =
        P->data[j * P->n + j];

        if (!isfinite(diagonal) ||
            diagonal <= 0.0)
        {
            return -1;
        }

        charge[j] =
        V_c[j] / diagonal;
    }

    return 0;
}


/*
 * ============================================================
 * phi = S V_c = P F V_c
 * ============================================================
 */
int potential_self_voltage_to_full_potential(
    const PotentialMatrix *P,
    const double *V_c,
    double *phi
)
{
    if (P == NULL ||
        P->data == NULL ||
        P->n == 0 ||
        V_c == NULL ||
        phi == NULL)
    {
        return -1;
    }

    size_t nv =
    P->n;

    for (size_t j = 0; j < nv; ++j)
    {
        double value =
        0.0;

        for (size_t a = 0; a < nv; ++a)
        {
            double P_aa =
            P->data[a * nv + a];

            if (!isfinite(P_aa) ||
                P_aa <= 0.0)
            {
                return -1;
            }

            value +=
            P->data[j * nv + a]
            * V_c[a]
            / P_aa;
        }

        phi[j] =
        value;
    }

    return 0;
}
