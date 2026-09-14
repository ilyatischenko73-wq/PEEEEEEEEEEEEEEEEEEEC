#ifndef PEEC_APERTURE_COUPLING_H
#define PEEC_APERTURE_COUPLING_H

#include "mesh.h"

#include <stddef.h>


/*
 * ============================================================
 * КЛАСС РЕБРА ВРЕМЕННОЙ PEC-КРЫШКИ
 * ============================================================
 *
 * PERIMETER:
 *
 *     оба конца лежат на границе раскрыва.
 *
 *
 * CUT:
 *
 *     один конец лежит на границе раскрыва,
 *     второй лежит внутри удаляемой PEC-крышки.
 *
 *
 * INTERNAL:
 *
 *     оба конца лежат внутри крышки.
 *
 *
 * Для Stage 1 -> Stage 2 источника используются именно
 * CUT edges.
 */
typedef enum
{
    APERTURE_EDGE_PERIMETER = 0,
    APERTURE_EDGE_CUT       = 1,
    APERTURE_EDGE_INTERNAL  = 2

} ApertureEdgeClass;


/*
 * ============================================================
 * СООТВЕТСТВИЕ УЗЛОВ CLOSED -> OPEN
 * ============================================================
 *
 * Один элемент задаёт:
 *
 *     closed boundary node
 *
 *         ->
 *
 *     open boundary node.
 *
 *
 * Именно это соответствие в Python читалось из
 * aperture_map.json.
 */
typedef struct
{
    size_t closed_node;
    size_t open_node;

} ApertureNodeMapEntry;


/*
 * ============================================================
 * ОПЕРАТОР СВЯЗИ STAGE 1 -> STAGE 2
 * ============================================================
 */
typedef struct
{
    /*
     * --------------------------------------------------------
     * РАЗМЕРЫ
     * --------------------------------------------------------
     */

    /*
     * Количество узлов открытой Stage-2 сетки.
     */
    size_t n_open_nodes;

    /*
     * Количество ребер временной PEC-крышки.
     */
    size_t n_cover_edges;

    /*
     * Количество correspondence entries:
     *
     *     closed boundary node -> open boundary node.
     */
    size_t n_node_map;


    /*
     * --------------------------------------------------------
     * NODE MAPPING
     * --------------------------------------------------------
     */
    ApertureNodeMapEntry *node_map;


    /*
     * --------------------------------------------------------
     * COVER EDGE TOPOLOGY
     * --------------------------------------------------------
     *
     * Для каждого локального ребра крышки:
     *
     *     edge_nodes[2*e + 0] = a
     *     edge_nodes[2*e + 1] = b
     *
     * где a,b — номера узлов CLOSED mesh.
     *
     *
     * Положительное направление тока:
     *
     *     a -> b.
     */
    size_t *edge_nodes;


    /*
     * --------------------------------------------------------
     * КЛАССИФИКАЦИЯ РЕБЕР
     * --------------------------------------------------------
     */
    ApertureEdgeClass *edge_class;

    size_t n_perimeter_edges;
    size_t n_cut_edges;
    size_t n_internal_edges;


    /*
     * --------------------------------------------------------
     * B_GAMMA
     * --------------------------------------------------------
     *
     * Матрица:
     *
     *     B_gamma
     *
     * размер:
     *
     *     Nv_open x Ne_cover.
     *
     *
     * Связь:
     *
     *     I_gamma = B_gamma I_cover.
     */
    double *B_gamma;


    /*
     * --------------------------------------------------------
     * ОПТИМИЗАЦИЯ ДЛЯ STAGE 2
     * --------------------------------------------------------
     *
     * Если P_open известна:
     *
     *     P_B_gamma = P_open B_gamma.
     *
     *
     * Тогда на каждом временном шаге:
     *
     *     P I_gamma
     *
     *     =
     *
     *     P_B_gamma I_cover.
     */
    double *P_B_gamma;

    int has_P_B_gamma;

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
 * СОЗДАНИЕ ИЗ TOPOLOGY + NODE MAP
 * ============================================================
 *
 * cover_edge_nodes:
 *
 *     массив длины:
 *
 *         2 * n_cover_edges.
 *
 *
 * node_map:
 *
 *     correspondence:
 *
 *         closed boundary node -> open boundary node.
 *
 *
 * Функция:
 *
 * 1. копирует topology крышки;
 * 2. копирует node map;
 * 3. классифицирует ребра;
 * 4. строит B_gamma.
 */
int aperture_coupling_build(
    ApertureCoupling *coupling,
    size_t n_open_nodes,
    size_t n_cover_edges,
    const size_t *cover_edge_nodes,
    size_t n_node_map,
    const ApertureNodeMapEntry *node_map
);


/*
 * ============================================================
 * ПОИСК CLOSED NODE В NODE MAP
 * ============================================================
 *
 * Возвращает:
 *
 *     1  если узел найден;
 *     0  если узел не найден;
 *    -1  при ошибке.
 */
int aperture_coupling_find_open_node(
    const ApertureCoupling *coupling,
    size_t closed_node,
    size_t *open_node
);


/*
 * ============================================================
 * КЛАССИФИКАЦИЯ РЕБЕР
 * ============================================================
 */
int aperture_coupling_classify_edges(
    ApertureCoupling *coupling
);


/*
 * ============================================================
 * BUILD B_GAMMA
 * ============================================================
 *
 * Для CUT edge:
 *
 *     a -> b.
 *
 *
 * Если a лежит на boundary:
 *
 *     B_gamma[open(a), e] = -1.
 *
 *
 * Если b лежит на boundary:
 *
 *     B_gamma[open(b), e] = +1.
 *
 *
 * Это полностью повторяет дискретную Stage-1 -> Stage-2
 * схему старой Python-версии.
 */
int aperture_coupling_build_B_gamma(
    ApertureCoupling *coupling
);


/*
 * ============================================================
 * I_GAMMA
 * ============================================================
 *
 * Вычисляет:
 *
 *     I_gamma =
 *
 *     B_gamma I_cover.
 *
 *
 * Размеры:
 *
 *     I_cover : Ne_cover
 *     I_gamma : Nv_open
 */
int aperture_coupling_compute_node_source(
    const ApertureCoupling *coupling,
    const double *cover_current,
    double *node_source
);


/*
 * ============================================================
 * NET CURRENT
 * ============================================================
 *
 * Вычисляет:
 *
 *     sum_j I_gamma[j].
 *
 *
 * Эта величина полезна как диагностическая проверка.
 */
double aperture_coupling_compute_net_current(
    const ApertureCoupling *coupling,
    const double *cover_current
);


/*
 * ============================================================
 * P B_GAMMA
 * ============================================================
 *
 * P_open:
 *
 *     Nv_open x Nv_open.
 *
 *
 * Формирует:
 *
 *     P_B_gamma =
 *
 *     P_open B_gamma.
 *
 *
 * Размер:
 *
 *     Nv_open x Ne_cover.
 */
int aperture_coupling_build_P_B_gamma(
    ApertureCoupling *coupling,
    const double *P_open
);


/*
 * ============================================================
 * P I_GAMMA
 * ============================================================
 *
 * Быстро вычисляет:
 *
 *     P I_gamma
 *
 *     =
 *
 *     P_B_gamma I_cover.
 *
 *
 * Требует предварительного вызова:
 *
 *     aperture_coupling_build_P_B_gamma().
 */
int aperture_coupling_compute_P_node_source(
    const ApertureCoupling *coupling,
    const double *cover_current,
    double *P_node_source
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
