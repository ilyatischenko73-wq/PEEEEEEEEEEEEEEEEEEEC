#ifndef PEEC_DUAL_MESH_H
#define PEEC_DUAL_MESH_H

#include "mesh.h"

#include <stddef.h>


/*
 * Один четырехугольный участок дуальной области.
 *
 * points:
 *     4 вершины quadrangle, каждая по 3 координаты.
 *
 * points[3*k + 0] = x_k
 * points[3*k + 1] = y_k
 * points[3*k + 2] = z_k
 *
 * owner_id:
 *     индекс узла или ребра, которому принадлежит патч.
 *
 * cell_id:
 *     индекс исходной четырехугольной ячейки.
 */
typedef struct
{
    double points[12];

    int owner_id;
    int cell_id;

} QuadPatch;


/*
 * Динамический список QuadPatch.
 *
 * Для каждого узла и каждого ребра число патчей различается:
 *
 * - внутренний узел обычно получает несколько quarter-quads;
 * - внутреннее ребро обычно получает два half-quads;
 * - граничное ребро получает один half-quad.
 */
typedef struct
{
    QuadPatch *patches;

    size_t count;
    size_t capacity;

} PatchList;


/*
 * Дуальная геометрия метода ЧЭС.
 *
 * node_patches[j]:
 *     геометрическая область Pi_j^v вокруг узла j.
 *
 * edge_patches[i]:
 *     геометрическая область Pi_i^e, связанная с ребром i.
 *
 * node_region_areas[j]:
 *     площадь Pi_j^v.
 *
 * edge_region_areas[i]:
 *     площадь Pi_i^e.
 *
 * node_region_centers:
 *     геометрические центры узловых областей.
 *
 * edge_region_centers:
 *     геометрические центры реберных областей.
 */
typedef struct
{
    size_t n_nodes;
    size_t n_edges;

    PatchList *node_patches;
    PatchList *edge_patches;

    double *node_region_areas;
    double *edge_region_areas;

    double *node_region_centers;
    double *edge_region_centers;

} DualMesh;


/*
 * Инициализация пустой структуры.
 */
void dual_mesh_init(DualMesh *dual);


/*
 * Освобождение всей памяти дуальной сетки.
 */
void dual_mesh_free(DualMesh *dual);


/*
 * Построение дуальной геометрии из исходной quad-сетки.
 *
 * Используется та же конструкция, что и в Python-версии:
 *
 * - узловая область строится из quarter-quads;
 * - реберная область строится из half-quads.
 */
int dual_mesh_build(const Mesh *mesh, DualMesh *dual);


/*
 * Печатает краткую диагностическую информацию.
 */
void dual_mesh_print_info(const DualMesh *dual);


#endif
