#include "dual_mesh.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * Норма трехмерного вектора.
 */
static double norm3(double x, double y, double z)
{
    return sqrt(x * x + y * y + z * z);
}


/*
 * Векторное произведение:
 *
 * result = a x b.
 */
static void cross3(const double a[3], const double b[3], double result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}


/*
 * Копирование одной 3D-точки.

static void copy_point(double dst[3], const double src[3])
{
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
}
!времменно не используется

*/
/*
 * Вычисляет середину отрезка:
 *
 * m = (a + b) / 2.
 */
static void midpoint(const double a[3], const double b[3], double m[3])
{
    m[0] = 0.5 * (a[0] + b[0]);
    m[1] = 0.5 * (a[1] + b[1]);
    m[2] = 0.5 * (a[2] + b[2]);
}


/*
 * Геометрический центр четырехугольника:
 *
 * c = (p0 + p1 + p2 + p3) / 4.
 */
static void quad_center(double p[4][3], double center[3])
{
    center[0] = 0.25 * (p[0][0] + p[1][0] + p[2][0] + p[3][0]);
    center[1] = 0.25 * (p[0][1] + p[1][1] + p[2][1] + p[3][1]);
    center[2] = 0.25 * (p[0][2] + p[1][2] + p[2][2] + p[3][2]);
}


/*
 * Площадь произвольного quadrangle.
 *
 * Разбиваем его на два треугольника:
 *
 * (p0, p1, p2)
 * (p0, p2, p3)
 *
 * Такой же подход использован для исходной quad-сетки.
 */
static double quad_area(double p[4][3])
{
    double v01[3] = {
        p[1][0] - p[0][0],
        p[1][1] - p[0][1],
        p[1][2] - p[0][2]
    };

    double v02[3] = {
        p[2][0] - p[0][0],
        p[2][1] - p[0][1],
        p[2][2] - p[0][2]
    };

    double v03[3] = {
        p[3][0] - p[0][0],
        p[3][1] - p[0][1],
        p[3][2] - p[0][2]
    };

    double c1[3];
    double c2[3];

    cross3(v01, v02, c1);
    cross3(v02, v03, c2);

    double a1 = 0.5 * norm3(c1[0], c1[1], c1[2]);
    double a2 = 0.5 * norm3(c2[0], c2[1], c2[2]);

    return a1 + a2;
}


/*
 * Центр quadrangle для взвешивания геометрического центра
 * составной области.
 */
static void patch_center(const QuadPatch *patch, double center[3])
{
    center[0] = 0.25 * (
        patch->points[0] +
        patch->points[3] +
        patch->points[6] +
        patch->points[9]
    );

    center[1] = 0.25 * (
        patch->points[1] +
        patch->points[4] +
        patch->points[7] +
        patch->points[10]
    );

    center[2] = 0.25 * (
        patch->points[2] +
        patch->points[5] +
        patch->points[8] +
        patch->points[11]
    );
}


/*
 * Переводит QuadPatch в обычный массив p[4][3].
 */
static void patch_to_points(const QuadPatch *patch, double p[4][3])
{
    for (int k = 0; k < 4; ++k)
    {
        p[k][0] = patch->points[3 * k + 0];
        p[k][1] = patch->points[3 * k + 1];
        p[k][2] = patch->points[3 * k + 2];
    }
}


/*
 * Добавление нового патча в динамический список.
 */
static int patch_list_append(PatchList *list, const QuadPatch *patch)
{
    if (list->count >= list->capacity)
    {
        size_t new_capacity = (list->capacity == 0) ? 4 : 2 * list->capacity;

        QuadPatch *new_patches = realloc(
            list->patches,
            new_capacity * sizeof(QuadPatch)
        );

        if (new_patches == NULL)
        {
            return -1;
        }

        list->patches = new_patches;
        list->capacity = new_capacity;
    }

    list->patches[list->count] = *patch;
    ++list->count;

    return 0;
}


/*
 * Заполняет один QuadPatch четырьмя точками.
 */
static void make_patch(
    QuadPatch *patch,
    const double p0[3],
    const double p1[3],
    const double p2[3],
    const double p3[3],
    int owner_id,
    int cell_id
)
{
    const double *points[4] = {p0, p1, p2, p3};

    for (int k = 0; k < 4; ++k)
    {
        patch->points[3 * k + 0] = points[k][0];
        patch->points[3 * k + 1] = points[k][1];
        patch->points[3 * k + 2] = points[k][2];
    }

    patch->owner_id = owner_id;
    patch->cell_id = cell_id;
}


void dual_mesh_init(DualMesh *dual)
{
    if (dual == NULL)
    {
        return;
    }

    memset(dual, 0, sizeof(*dual));
}


void dual_mesh_free(DualMesh *dual)
{
    if (dual == NULL)
    {
        return;
    }

    /*
     * Освобождаем все списки узловых патчей.
     */
    if (dual->node_patches != NULL)
    {
        for (size_t i = 0; i < dual->n_nodes; ++i)
        {
            free(dual->node_patches[i].patches);
        }
    }

    /*
     * Освобождаем все списки реберных патчей.
     */
    if (dual->edge_patches != NULL)
    {
        for (size_t i = 0; i < dual->n_edges; ++i)
        {
            free(dual->edge_patches[i].patches);
        }
    }

    free(dual->node_patches);
    free(dual->edge_patches);

    free(dual->node_region_areas);
    free(dual->edge_region_areas);

    free(dual->node_region_centers);
    free(dual->edge_region_centers);

    dual_mesh_init(dual);
}


/*
 * Добавляет quarter-quad в область конкретного узла.
 *
 * Для исходной ячейки:
 *
 *      q3 -------- q2
 *       |          |
 *       |          |
 *      q0 -------- q1
 *
 * строятся середины ребер:
 *
 * m01, m12, m23, m30
 *
 * и центр c.
 *
 * Узловая область q0:
 *
 *      m30 ----- c
 *       |        |
 *       |        |
 *      q0 ----- m01
 *
 * то есть:
 *
 * [q0, m01, c, m30]
 */
static int add_node_patch(
    DualMesh *dual,
    int node_id,
    int cell_id,
    const double p0[3],
    const double p1[3],
    const double p2[3],
    const double p3[3]
)
{
    QuadPatch patch;

    make_patch(
        &patch,
        p0,
        p1,
        p2,
        p3,
        node_id,
        cell_id
    );

    return patch_list_append(
        &dual->node_patches[node_id],
        &patch
    );
}


/*
 * Добавляет half-quad в область конкретного ребра.
 *
 * В нашей ЧЭС-нормировке область Pi_i^e строится отдельно
 * в каждой соседней ячейке.
 */
static int add_edge_patch(
    DualMesh *dual,
    int edge_id,
    int cell_id,
    const double p0[3],
    const double p1[3],
    const double p2[3],
    const double p3[3]
)
{
    QuadPatch patch;

    make_patch(
        &patch,
        p0,
        p1,
        p2,
        p3,
        edge_id,
        cell_id
    );

    return patch_list_append(
        &dual->edge_patches[edge_id],
        &patch
    );
}


int dual_mesh_build(const Mesh *mesh, DualMesh *dual)
{
    if (mesh == NULL || dual == NULL)
    {
        return -1;
    }

    if (mesh->xyz == NULL ||
        mesh->quads == NULL ||
        mesh->cell_edges == NULL)
    {
        fprintf(stderr, "ERROR: topology must be built before dual mesh.\n");
        return -1;
    }

    dual_mesh_init(dual);

    dual->n_nodes = mesh->n_nodes;
    dual->n_edges = mesh->n_edges;

    /*
     * Для каждого узла и ребра создаем отдельный динамический список патчей.
     */
    dual->node_patches = calloc(mesh->n_nodes, sizeof(PatchList));
    dual->edge_patches = calloc(mesh->n_edges, sizeof(PatchList));

    dual->node_region_areas = calloc(mesh->n_nodes, sizeof(double));
    dual->edge_region_areas = calloc(mesh->n_edges, sizeof(double));

    dual->node_region_centers = calloc(3 * mesh->n_nodes, sizeof(double));
    dual->edge_region_centers = calloc(3 * mesh->n_edges, sizeof(double));

    if (dual->node_patches == NULL ||
        dual->edge_patches == NULL ||
        dual->node_region_areas == NULL ||
        dual->edge_region_areas == NULL ||
        dual->node_region_centers == NULL ||
        dual->edge_region_centers == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate dual mesh.\n");
        dual_mesh_free(dual);
        return -1;
    }

    /*
     * =======================================================
     * ПРОХОД ПО ВСЕМ ИСХОДНЫМ QUAD-ЯЧЕЙКАМ
     * =======================================================
     */
    for (size_t cell = 0; cell < mesh->n_quads; ++cell)
    {
        const int *q = &mesh->quads[4 * cell];

        /*
         * Координаты четырех вершин исходной ячейки.
         */
        double p[4][3];

        for (int k = 0; k < 4; ++k)
        {
            p[k][0] = mesh->xyz[3 * q[k] + 0];
            p[k][1] = mesh->xyz[3 * q[k] + 1];
            p[k][2] = mesh->xyz[3 * q[k] + 2];
        }

        /*
         * Середины четырех ребер ячейки.
         */
        double m01[3];
        double m12[3];
        double m23[3];
        double m30[3];

        midpoint(p[0], p[1], m01);
        midpoint(p[1], p[2], m12);
        midpoint(p[2], p[3], m23);
        midpoint(p[3], p[0], m30);

        /*
         * Геометрический центр исходной quad-ячейки.
         */
        double center[3];

        quad_center(p, center);

        /*
         * ===================================================
         * УЗЛОВЫЕ ОБЛАСТИ Pi_j^v
         * ===================================================
         *
         * В каждой исходной ячейке формируются четыре quarter-quads:
         *
         * q0 : [p0, m01, center, m30]
         * q1 : [p1, m12, center, m01]
         * q2 : [p2, m23, center, m12]
         * q3 : [p3, m30, center, m23]
         */

        if (add_node_patch(
            dual,
            q[0],
            (int)cell,
                           p[0],
                           m01,
                           center,
                           m30
        ) != 0)
        {
            goto allocation_error;
        }

        if (add_node_patch(
            dual,
            q[1],
            (int)cell,
                           p[1],
                           m12,
                           center,
                           m01
        ) != 0)
        {
            goto allocation_error;
        }

        if (add_node_patch(
            dual,
            q[2],
            (int)cell,
                           p[2],
                           m23,
                           center,
                           m12
        ) != 0)
        {
            goto allocation_error;
        }

        if (add_node_patch(
            dual,
            q[3],
            (int)cell,
                           p[3],
                           m30,
                           center,
                           m23
        ) != 0)
        {
            goto allocation_error;
        }

        /*
         * ===================================================
         * РЕБЕРНЫЕ ОБЛАСТИ Pi_i^e
         * ===================================================
         *
         * Для каждой стороны исходной quad-ячейки добавляем
         * один half-quad.
         *
         * Схема совпадает с нашей Python-реализацией:
         *
         * edge 0, q0-q1:
         *
         *     [p0, p1, m12, m30]
         *
         * edge 1, q1-q2:
         *
         *     [p1, p2, m23, m01]
         *
         * edge 2, q2-q3:
         *
         *     [p2, p3, m30, m12]
         *
         * edge 3, q3-q0:
         *
         *     [p3, p0, m01, m23]
         */

        int e0 = mesh->cell_edges[4 * cell + 0];
        int e1 = mesh->cell_edges[4 * cell + 1];
        int e2 = mesh->cell_edges[4 * cell + 2];
        int e3 = mesh->cell_edges[4 * cell + 3];

        if (add_edge_patch(
            dual,
            e0,
            (int)cell,
                           p[0],
                           p[1],
                           m12,
                           m30
        ) != 0)
        {
            goto allocation_error;
        }

        if (add_edge_patch(
            dual,
            e1,
            (int)cell,
                           p[1],
                           p[2],
                           m23,
                           m01
        ) != 0)
        {
            goto allocation_error;
        }

        if (add_edge_patch(
            dual,
            e2,
            (int)cell,
                           p[2],
                           p[3],
                           m30,
                           m12
        ) != 0)
        {
            goto allocation_error;
        }

        if (add_edge_patch(
            dual,
            e3,
            (int)cell,
                           p[3],
                           p[0],
                           m01,
                           m23
        ) != 0)
        {
            goto allocation_error;
        }
    }

    /*
     * =======================================================
     * ПЛОЩАДИ И ЦЕНТРЫ УЗЛОВЫХ ОБЛАСТЕЙ
     * =======================================================
     */
    for (size_t node = 0; node < dual->n_nodes; ++node)
    {
        PatchList *list = &dual->node_patches[node];

        double total_area = 0.0;
        double weighted_center[3] = {0.0, 0.0, 0.0};

        for (size_t k = 0; k < list->count; ++k)
        {
            QuadPatch *patch = &list->patches[k];

            double p[4][3];
            double center[3];

            patch_to_points(patch, p);
            patch_center(patch, center);

            double area = quad_area(p);

            total_area += area;

            weighted_center[0] += area * center[0];
            weighted_center[1] += area * center[1];
            weighted_center[2] += area * center[2];
        }

        dual->node_region_areas[node] = total_area;

        if (total_area > 0.0)
        {
            dual->node_region_centers[3 * node + 0] =
            weighted_center[0] / total_area;

            dual->node_region_centers[3 * node + 1] =
            weighted_center[1] / total_area;

            dual->node_region_centers[3 * node + 2] =
            weighted_center[2] / total_area;
        }
    }

    /*
     * =======================================================
     * ПЛОЩАДИ И ЦЕНТРЫ РЕБЕРНЫХ ОБЛАСТЕЙ
     * =======================================================
     */
    for (size_t edge = 0; edge < dual->n_edges; ++edge)
    {
        PatchList *list = &dual->edge_patches[edge];

        double total_area = 0.0;
        double weighted_center[3] = {0.0, 0.0, 0.0};

        for (size_t k = 0; k < list->count; ++k)
        {
            QuadPatch *patch = &list->patches[k];

            double p[4][3];
            double center[3];

            patch_to_points(patch, p);
            patch_center(patch, center);

            double area = quad_area(p);

            total_area += area;

            weighted_center[0] += area * center[0];
            weighted_center[1] += area * center[1];
            weighted_center[2] += area * center[2];
        }

        dual->edge_region_areas[edge] = total_area;

        if (total_area > 0.0)
        {
            dual->edge_region_centers[3 * edge + 0] =
            weighted_center[0] / total_area;

            dual->edge_region_centers[3 * edge + 1] =
            weighted_center[1] / total_area;

            dual->edge_region_centers[3 * edge + 2] =
            weighted_center[2] / total_area;
        }
    }

    return 0;


    allocation_error:

    fprintf(stderr, "ERROR: cannot allocate dual mesh patch.\n");
    dual_mesh_free(dual);

    return -1;
}


void dual_mesh_print_info(const DualMesh *dual)
{
    if (dual == NULL)
    {
        return;
    }

    size_t total_node_patches = 0;
    size_t total_edge_patches = 0;

    size_t min_node_patches = (size_t)-1;
    size_t max_node_patches = 0;

    size_t min_edge_patches = (size_t)-1;
    size_t max_edge_patches = 0;

    for (size_t i = 0; i < dual->n_nodes; ++i)
    {
        size_t count = dual->node_patches[i].count;

        total_node_patches += count;

        if (count < min_node_patches)
        {
            min_node_patches = count;
        }

        if (count > max_node_patches)
        {
            max_node_patches = count;
        }
    }

    for (size_t i = 0; i < dual->n_edges; ++i)
    {
        size_t count = dual->edge_patches[i].count;

        total_edge_patches += count;

        if (count < min_edge_patches)
        {
            min_edge_patches = count;
        }

        if (count > max_edge_patches)
        {
            max_edge_patches = count;
        }
    }

    printf("\n");
    printf("========================================\n");
    printf("DUAL MESH\n");
    printf("========================================\n");
    printf("Node regions       : %zu\n", dual->n_nodes);
    printf("Edge regions       : %zu\n", dual->n_edges);
    printf("Node patches total : %zu\n", total_node_patches);
    printf("Edge patches total : %zu\n", total_edge_patches);
    printf("Node patches min   : %zu\n", min_node_patches);
    printf("Node patches max   : %zu\n", max_node_patches);
    printf("Edge patches min   : %zu\n", min_edge_patches);
    printf("Edge patches max   : %zu\n", max_edge_patches);
    printf("========================================\n");
}
