#include "inductance.h"
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
 * Магнитная постоянная:
 *
 *     mu0 = 4*pi*10^-7 Гн/м.
 *
 * Поэтому:
 *
 *     mu0/(4*pi) = 10^-7 Гн/м.
 */
static const double MU0_OVER_4PI = 1.0e-7;


/*
 * ============================================================
 * СКАЛЯРНОЕ ПРОИЗВЕДЕНИЕ
 * ============================================================
 */
static double dot3(const double a[3], const double b[3])
{
    return
    a[0] * b[0]
    + a[1] * b[1]
    + a[2] * b[2];
}


/*
 * ============================================================
 * ПАРАМЕТРЫ ПО УМОЛЧАНИЮ
 * ============================================================
 */
void inductance_options_default(
    InductanceOptions *options,
    int parallel_threads
)
{
    if (options == NULL)
    {
        return;
    }

    /*
     * Хорошо разделенные области.
     */
    options->far_order = 4;

    /*
     * Близкие, но несоприкасающиеся области.
     */
    options->near_order = 10;

    /*
     * Области с общей вершиной или ребром.
     *
     * Пока используется обычная квадратура повышенного порядка.
     */
    options->touching_order = 12;

    /*
     * SELF-интеграл.
     *
     * Сингулярность устраняется Duffy-преобразованием,
     * поэтому обычно достаточно умеренного порядка.
     */
    options->self_order = 8;

    /*
     * Критерий NEAR/FAR.
     */
    options->near_factor = 0.5;

    /*
     * Число OpenMP-потоков.
     */
    options->parallel_threads = parallel_threads;

    /*
     * Допуск при проверке:
     *
     *     e_a . e_b = 0.
     */
    options->direction_tolerance = 1.0e-10;
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void inductance_init(InductanceMatrix *L)
{
    if (L == NULL)
    {
        return;
    }

    memset(L, 0, sizeof(*L));
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ ПАМЯТИ
 * ============================================================
 */
void inductance_free(InductanceMatrix *L)
{
    if (L == NULL)
    {
        return;
    }

    free(L->data);

    inductance_init(L);
}


/*
 * ============================================================
 * ЭФФЕКТИВНАЯ ШИРИНА РЕБЕРНОЙ ОБЛАСТИ
 * ============================================================
 *
 * В нашей нормировке:
 *
 *     J_i(r,t) = I_i(t)/w_i * e_i.
 *
 * При этом:
 *
 *     w_i = S_i^e / l_i.
 *
 * Здесь:
 *
 *     S_i^e - площадь Pi_i^e;
 *     l_i   - длина глобального ребра.
 */
static double edge_effective_width(
    const Mesh *mesh,
    const DualMesh *dual,
    size_t edge
)
{
    double area = dual->edge_region_areas[edge];
    double length = mesh->edge_lengths[edge];

    if (area <= 0.0 || length <= 0.0)
    {
        return 0.0;
    }

    return area / length;
}


/*
 * ============================================================
 * ИНТЕГРАЛ ДЛЯ ДВУХ QUADPATCH
 * ============================================================
 *
 * quadrature.c самостоятельно выбирает:
 *
 *     SELF
 *     TOUCHING
 *     NEAR
 *     FAR.
 */
static int compute_patch_pair_kernel(
    const QuadPatch *patch_a,
    const QuadPatch *patch_b,
    const InductanceOptions *options,
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
 * ОДИН КОЭФФИЦИЕНТ L_ab
 * ============================================================
 */
int inductance_compute_pair(
    const Mesh *mesh,
    const DualMesh *dual,
    size_t edge_a,
    size_t edge_b,
    const InductanceOptions *options,
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

    if (edge_a >= mesh->n_edges || edge_b >= mesh->n_edges)
    {
        return -1;
    }

    /*
     * Единичные направления глобальных ребер:
     *
     *     e_a
     *     e_b.
     */
    const double *direction_a =
    &mesh->edge_directions[3 * edge_a];

    const double *direction_b =
    &mesh->edge_directions[3 * edge_b];

    /*
     * Ориентационный множитель:
     *
     *     e_a . e_b.
     */
    double direction_dot = dot3(
        direction_a,
        direction_b
    );

    /*
     * Если ребра ортогональны, в принятой базисной модели:
     *
     *     L_ab = 0.
     *
     * Эту проверку выполняем до дорогого интегрирования.
     */
    if (fabs(direction_dot) < options->direction_tolerance)
    {
        *value = 0.0;

        return 0;
    }

    /*
     * Эффективные ширины:
     *
     *     w_a = S_a^e/l_a
     *     w_b = S_b^e/l_b.
     */
    double width_a = edge_effective_width(
        mesh,
        dual,
        edge_a
    );

    double width_b = edge_effective_width(
        mesh,
        dual,
        edge_b
    );

    if (width_a <= 0.0 || width_b <= 0.0)
    {
        return -1;
    }

    /*
     * Полные реберные области Pi_a^e и Pi_b^e
     * состоят из нескольких QuadPatch.
     */
    const PatchList *region_a =
    &dual->edge_patches[edge_a];

    const PatchList *region_b =
    &dual->edge_patches[edge_b];

    if (region_a->count == 0 || region_b->count == 0)
    {
        return -1;
    }

    /*
     * Полный геометрический интеграл:
     *
     *                  /        /
     *                 |        |       1
     *     G_ab =      |        |    ------- dS dS'.
     *                 |        |    |r-r'|
     *                /Pi_a    /Pi_b
     *
     * Поскольку Pi_a и Pi_b составные,
     * суммируем все пары QuadPatch.
     */
    double integral = 0.0;

    for (size_t i = 0; i < region_a->count; ++i)
    {
        const QuadPatch *patch_a =
        &region_a->patches[i];

        for (size_t j = 0; j < region_b->count; ++j)
        {
            const QuadPatch *patch_b =
            &region_b->patches[j];

            double patch_value = 0.0;

            int status = compute_patch_pair_kernel(
                patch_a,
                patch_b,
                options,
                &patch_value
            );

            if (status != 0)
            {
                return -1;
            }

            integral += patch_value;
        }
    }

    /*
     * Полная формула:
     *
     *              mu0    e_a . e_b
     *     L_ab = -------  ----------- G_ab.
     *              4*pi      w_a w_b
     *
     * Размерность:
     *
     *     Гн.
     */
    *value =
    MU0_OVER_4PI
    * direction_dot
    * integral
    / (width_a * width_b);

    if (!isfinite(*value))
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ПОЛНАЯ МАТРИЦА L
 * ============================================================
 */
int inductance_compute_matrix(
    const Mesh *mesh,
    const DualMesh *dual,
    const InductanceOptions *options,
    InductanceMatrix *L
)
{
    if (mesh == NULL ||
        dual == NULL ||
        options == NULL ||
        L == NULL)
    {
        return -1;
    }

    if (mesh->n_edges == 0)
    {
        return -1;
    }

    inductance_init(L);

    size_t ne = mesh->n_edges;

    L->n = ne;

    /*
     * Матрица плотная:
     *
     *     Ne x Ne.
     *
     * Используем calloc(), чтобы начальное содержимое
     * было гарантированно нулевым.
     */
    L->data = calloc(
        ne * ne,
        sizeof(double)
    );

    if (L->data == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate inductance matrix.\n"
        );

        return -1;
    }

    /*
     * Общий код ошибки для OpenMP region.
     *
     * Нельзя выполнять return непосредственно
     * из parallel for.
     */
    int error_code = 0;


    /*
     * ========================================================
     * OPENMP
     * ========================================================
     *
     * Матрица L симметрична:
     *
     *     L_ab = L_ba.
     *
     * Поэтому вычисляем только:
     *
     *     b >= a.
     *
     *
     * schedule(dynamic) полезен, потому что стоимость
     * различных пар неодинакова:
     *
     *     FAR       - дешевле;
     *     NEAR      - дороже;
     *     TOUCHING  - еще дороже;
     *     SELF      - отдельная Duffy-квадратура.
     */
    #pragma omp parallel for schedule(dynamic) num_threads(options->parallel_threads)
    for (long long a_signed = 0; a_signed < (long long)ne; ++a_signed)
    {
        size_t a = (size_t)a_signed;

        /*
         * Проверяем, не обнаружил ли другой поток ошибку.
         */
        int current_error = 0;

        #pragma omp atomic read
        current_error = error_code;

        if (current_error != 0)
        {
            continue;
        }

        for (size_t b = a; b < ne; ++b)
        {
            /*
             * Повторно проверяем флаг ошибки,
             * чтобы длинная строка могла завершиться раньше.
             */
            #pragma omp atomic read
            current_error = error_code;

            if (current_error != 0)
            {
                break;
            }

            double value = 0.0;

            int status = inductance_compute_pair(
                mesh,
                dual,
                a,
                b,
                options,
                &value
            );

            if (status != 0)
            {
                /*
                 * Ошибочная ветка встречается редко,
                 * поэтому critical практически не влияет
                 * на производительность.
                 */
                #pragma omp critical(inductance_error)
                {
                    if (error_code == 0)
                    {
                        error_code = status;

                        fprintf(
                            stderr,
                            "\nERROR: inductance calculation failed "
                            "for edges %zu, %zu.\n",
                            a,
                            b
                        );
                    }
                }

                break;
            }

            /*
             * Одна пара (a,b) вычисляется только одним потоком.
             *
             * Поэтому запись в эти две ячейки безопасна:
             *
             *     L_ab
             *     L_ba.
             */
            L->data[a * ne + b] = value;
            L->data[b * ne + a] = value;
        }
    }

    if (error_code != 0)
    {
        inductance_free(L);

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ДИАГНОСТИКА МАТРИЦЫ
 * ============================================================
 */
void inductance_print_info(const InductanceMatrix *L)
{
    if (L == NULL || L->data == NULL || L->n == 0)
    {
        return;
    }

    double min_value = L->data[0];
    double max_value = L->data[0];

    double min_diagonal = HUGE_VAL;
    double max_diagonal = -HUGE_VAL;

    double max_symmetry_error = 0.0;

    /*
     * Проверяем значения матрицы и симметрию.
     */
    for (size_t i = 0; i < L->n; ++i)
    {
        double diagonal = L->data[i * L->n + i];

        if (diagonal < min_diagonal)
        {
            min_diagonal = diagonal;
        }

        if (diagonal > max_diagonal)
        {
            max_diagonal = diagonal;
        }

        for (size_t j = 0; j < L->n; ++j)
        {
            double value = L->data[i * L->n + j];

            if (value < min_value)
            {
                min_value = value;
            }

            if (value > max_value)
            {
                max_value = value;
            }

            double symmetry_error = fabs(
                value - L->data[j * L->n + i]
            );

            if (symmetry_error > max_symmetry_error)
            {
                max_symmetry_error = symmetry_error;
            }
        }
    }

    /*
     * Объем памяти плотной double-матрицы.
     */
    double memory_mib =
    ((double)L->n
    * (double)L->n
    * (double)sizeof(double))
    / (1024.0 * 1024.0);

    printf("\n");
    printf("========================================\n");
    printf("INDUCTANCE MATRIX\n");
    printf("========================================\n");
    printf("Size               : %zu x %zu\n", L->n, L->n);
    printf("Memory             : %.3f MiB\n", memory_mib);
    printf("Minimum L          : %.9e H\n", min_value);
    printf("Maximum L          : %.9e H\n", max_value);
    printf("Minimum diagonal   : %.9e H\n", min_diagonal);
    printf("Maximum diagonal   : %.9e H\n", max_diagonal);
    printf("Max symmetry error : %.3e H\n", max_symmetry_error);
    printf("========================================\n");
}
