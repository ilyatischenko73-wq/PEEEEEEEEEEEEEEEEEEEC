#ifndef PEEC_APERTURE_COUPLING_H
#define PEEC_APERTURE_COUPLING_H

#include "mesh.h"
#include "potential.h"

#include <stddef.h>


/*
 * ============================================================
 * NODE MAPPING
 * ============================================================
 *
 * Один узел границы временной PEC-крышки Stage 1
 * соответствует одному узлу границы апертуры Stage 2.
 */
typedef struct
{
    size_t closed_node;
    size_t open_node;

} ApertureNodeMapEntry;


/*
 * ============================================================
 * CUT EDGE
 * ============================================================
 *
 * Перенос Stage 1 -> Stage 2 выполняется только через
 * cut edges временной PEC-крышки.
 *
 *
 * Для каждого такого ребра заранее храним:
 *
 *     closed_edge
 *
 *         номер ребра на закрытой Stage-1 сетке;
 *
 *
 *     open_node
 *
 *         узел открытой Stage-2 сетки, куда переносится
 *         вклад этого ребра;
 *
 *
 *     sign
 *
 *         знак в B_gamma.
 *
 *
 * При глобальной ориентации ребра:
 *
 *     a -> b
 *
 * и incidence convention:
 *
 *     A[e,a] = -1
 *     A[e,b] = +1
 *
 * имеем:
 *
 *     если boundary endpoint = a:
 *
 *         sign = -1;
 *
 *     если boundary endpoint = b:
 *
 *         sign = +1.
 */
typedef struct
{
    size_t closed_edge;
    size_t open_node;
    double sign;

} ApertureCutEdge;


/*
 * ============================================================
 * APERTURE COUPLING
 * ============================================================
 *
 * Реализует дискретный перенос:
 *
 *     I_gamma = B_gamma I_cover.
 *
 *
 * Это инженерная дискретная аппроксимация переноса токов
 * временной PEC-крышки Stage 1 в эквивалентный узловой
 * источник Stage 2.
 *
 *
 * Cover-edge classification:
 *
 *     perimeter:
 *         оба конца ребра принадлежат aperture boundary;
 *
 *     internal:
 *         ни один конец не принадлежит aperture boundary;
 *
 *     cut:
 *         ровно один конец принадлежит aperture boundary.
 *
 *
 * В B_gamma участвуют только cut edges.
 */
typedef struct
{
    /*
     * Размеры исходных сеток.
     */
    size_t n_closed_nodes;
    size_t n_closed_edges;
    size_t n_open_nodes;

    /*
     * closed boundary node -> open boundary node.
     */
    ApertureNodeMapEntry *node_mapping;
    size_t n_node_mapping;

    /*
     * Рёбра временной PEC-крышки на Stage-1 сетке.
     */
    size_t *cover_edges;
    size_t n_cover_edges;

    /*
     * Предварительно построенное разреженное представление
     * B_gamma: одна запись на cut edge.
     */
    ApertureCutEdge *cut_edges;
    size_t n_cut_edges;

    /*
     * Диагностика классификации.
     */
    size_t n_perimeter_edges;
    size_t n_internal_edges;

    int built;
    int initialized;

} ApertureCoupling;


/*
 * ============================================================
 * INIT / FREE
 * ============================================================
 */
void aperture_coupling_init(
    ApertureCoupling *coupling
);

void aperture_coupling_free(
    ApertureCoupling *coupling
);


/*
 * ============================================================
 * BUILD FROM ARRAYS
 * ============================================================
 *
 * node_mapping:
 *
 *     массив соответствий closed_node -> open_node.
 *
 *
 * cover_edges:
 *
 *     номера рёбер временной PEC-крышки на closed_mesh.
 *
 *
 * После вызова функция:
 *
 *     1. проверяет mapping;
 *     2. классифицирует cover edges;
 *     3. строит sparse B_gamma для cut edges.
 */
int aperture_coupling_build(
    ApertureCoupling *coupling,
    const Mesh *closed_mesh,
    const Mesh *open_mesh,
    const ApertureNodeMapEntry *node_mapping,
    size_t n_node_mapping,
    const size_t *cover_edges,
    size_t n_cover_edges
);


/*
 * ============================================================
 * LOAD MAP FILE
 * ============================================================
 *
 * Поддерживаемый JSON-формат:
 *
 * {
 *   "node_mapping": [
 *     {"closed_node": 10, "open_node": 7},
 *     {"closed_node": 11, "open_node": 8}
 *   ],
 *
 *   "cover_edges": [
 *     100, 101, 102, 103
 *   ]
 * }
 *
 *
 * Парсер намеренно небольшой и читает только эти два поля.
 *
 * Если существующий aperture_map.json имеет другой формат,
 * нужно будет изменить только эту функцию; физика переноса
 * и структура B_gamma останутся без изменений.
 */
int aperture_coupling_load_json(
    ApertureCoupling *coupling,
    const Mesh *closed_mesh,
    const Mesh *open_mesh,
    const char *filename
);


/*
 * ============================================================
 * VALIDATE
 * ============================================================
 */
int aperture_coupling_validate(
    const ApertureCoupling *coupling
);


/*
 * ============================================================
 * I_gamma = B_gamma I_cover
 * ============================================================
 *
 * closed_edge_current:
 *
 *     полный вектор Stage-1 токов по closed_mesh:
 *
 *         Ne_closed.
 *
 *
 * open_node_current:
 *
 *     выходной эквивалентный узловой источник:
 *
 *         Nv_open.
 *
 *
 * Функция сама обнуляет open_node_current.
 */
int aperture_coupling_apply(
    const ApertureCoupling *coupling,
    const double *closed_edge_current,
    double *open_node_current
);


/*
 * ============================================================
 * P I_gamma
 * ============================================================
 *
 * Строит сразу потенциально-взвешенный источник:
 *
 *     P_open I_gamma.
 *
 *
 * Это можно использовать для оптимизированной Stage-2
 * реализации без полного умножения P на узловой источник
 * на каждом временном шаге.
 */
int aperture_coupling_apply_potential_source(
    const ApertureCoupling *coupling,
    const PotentialMatrix *P_open,
    const double *closed_edge_current,
    double *potential_source
);


/*
 * ============================================================
 * DIAGNOSTICS
 * ============================================================
 */
void aperture_coupling_print_info(
    const ApertureCoupling *coupling
);


#endif
