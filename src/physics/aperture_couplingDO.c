#include "aperture_coupling.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * ============================================================
 * DENSE MATRIX INDEX
 * ============================================================
 */
static size_t matrix_index(
    size_t row,
    size_t col,
    size_t n_cols
)
{
    return row * n_cols + col;
}


/*
 * ============================================================
 * INIT
 * ============================================================
 */
void aperture_coupling_init(
    ApertureCoupling *coupling
)
{
    if (coupling == NULL)
    {
        return;
    }

    memset(
        coupling,
        0,
        sizeof(*coupling)
    );
}


/*
 * ============================================================
 * FREE
 * ============================================================
 */
void aperture_coupling_free(
    ApertureCoupling *coupling
)
{
    if (coupling == NULL)
    {
        return;
    }

    free(
        coupling->node_map
    );

    free(
        coupling->edge_nodes
    );

    free(
        coupling->edge_class
    );

    free(
        coupling->B_gamma
    );

    free(
        coupling->P_B_gamma
    );

    aperture_coupling_init(
        coupling
    );
}


/*
 * ============================================================
 * CLOSED NODE -> OPEN NODE
 * ============================================================
 */
int aperture_coupling_find_open_node(
    const ApertureCoupling *coupling,
    size_t closed_node,
    size_t *open_node
)
{
    if (coupling == NULL ||
        open_node == NULL ||
        coupling->node_map == NULL)
    {
        return -1;
    }

    for (size_t i = 0; i < coupling->n_node_map; ++i)
    {
        if (coupling->node_map[i].closed_node == closed_node)
        {
            *open_node =
            coupling->node_map[i].open_node;

            return 1;
        }
    }

    return 0;
}


/*
 * ============================================================
 * CLASSIFY EDGES
 * ============================================================
 */
int aperture_coupling_classify_edges(
    ApertureCoupling *coupling
)
{
    if (coupling == NULL ||
        coupling->edge_nodes == NULL ||
        coupling->edge_class == NULL ||
        coupling->node_map == NULL)
    {
        return -1;
    }

    coupling->n_perimeter_edges =
    0;

    coupling->n_cut_edges =
    0;

    coupling->n_internal_edges =
    0;


    for (size_t e = 0; e < coupling->n_cover_edges; ++e)
    {
        size_t a =
        coupling->edge_nodes[2 * e + 0];

        size_t b =
        coupling->edge_nodes[2 * e + 1];


        size_t open_a =
        0;

        size_t open_b =
        0;


        int a_boundary =
        aperture_coupling_find_open_node(
            coupling,
            a,
            &open_a
        );

        int b_boundary =
        aperture_coupling_find_open_node(
            coupling,
            b,
            &open_b
        );


        if (a_boundary < 0 ||
            b_boundary < 0)
        {
            return -1;
        }


        /*
         * ----------------------------------------------------
         * PERIMETER
         * ----------------------------------------------------
         *
         * Оба узла ребра лежат на boundary.
         */
        if (a_boundary == 1 &&
            b_boundary == 1)
        {
            coupling->edge_class[e] =
            APERTURE_EDGE_PERIMETER;

            coupling->n_perimeter_edges++;

            continue;
        }


        /*
         * ----------------------------------------------------
         * INTERNAL
         * ----------------------------------------------------
         *
         * Ни один узел не лежит на boundary.
         */
        if (a_boundary == 0 &&
            b_boundary == 0)
        {
            coupling->edge_class[e] =
            APERTURE_EDGE_INTERNAL;

            coupling->n_internal_edges++;

            continue;
        }


        /*
         * ----------------------------------------------------
         * CUT
         * ----------------------------------------------------
         *
         * Ровно один узел лежит на boundary.
         */
        coupling->edge_class[e] =
        APERTURE_EDGE_CUT;

        coupling->n_cut_edges++;
    }


    return 0;
}


/*
 * ============================================================
 * BUILD B_GAMMA
 * ============================================================
 */
int aperture_coupling_build_B_gamma(
    ApertureCoupling *coupling
)
{
    if (coupling == NULL ||
        coupling->B_gamma == NULL ||
        coupling->edge_nodes == NULL ||
        coupling->edge_class == NULL)
    {
        return -1;
    }


    size_t nv =
    coupling->n_open_nodes;

    size_t ne =
    coupling->n_cover_edges;


    memset(
        coupling->B_gamma,
        0,
        nv * ne * sizeof(double)
    );


    for (size_t e = 0; e < ne; ++e)
    {
        if (coupling->edge_class[e] !=
            APERTURE_EDGE_CUT)
        {
            continue;
        }


        size_t a =
        coupling->edge_nodes[2 * e + 0];

        size_t b =
        coupling->edge_nodes[2 * e + 1];


        size_t open_a =
        0;

        size_t open_b =
        0;


        int a_boundary =
        aperture_coupling_find_open_node(
            coupling,
            a,
            &open_a
        );

        int b_boundary =
        aperture_coupling_find_open_node(
            coupling,
            b,
            &open_b
        );


        if (a_boundary < 0 ||
            b_boundary < 0)
        {
            return -1;
        }


        /*
         * ----------------------------------------------------
         * Положительный ток ребра:
         *
         *     a -> b.
         *
         *
         * Узловая инцидентность:
         *
         *     a : -I
         *     b : +I.
         * ----------------------------------------------------
         */


        /*
         * Boundary находится на конце a.
         */
        if (a_boundary == 1 &&
            b_boundary == 0)
        {
            if (open_a >= nv)
            {
                return -1;
            }

            coupling->B_gamma[
                matrix_index(
                    open_a,
                    e,
                    ne
                )
            ] =
            -1.0;

            continue;
        }


        /*
         * Boundary находится на конце b.
         */
        if (a_boundary == 0 &&
            b_boundary == 1)
        {
            if (open_b >= nv)
            {
                return -1;
            }

            coupling->B_gamma[
                matrix_index(
                    open_b,
                    e,
                    ne
                )
            ] =
            +1.0;

            continue;
        }


        /*
         * Если мы сюда попали, классификация и topology
         * противоречат друг другу.
         */
        fprintf(
            stderr,
            "ERROR: inconsistent CUT edge %zu.\n",
            e
        );

        return -1;
    }


    return 0;
}


/*
 * ============================================================
 * BUILD
 * ============================================================
 */
int aperture_coupling_build(
    ApertureCoupling *coupling,
    size_t n_open_nodes,
    size_t n_cover_edges,
    const size_t *cover_edge_nodes,
    size_t n_node_map,
    const ApertureNodeMapEntry *node_map
)
{
    if (coupling == NULL ||
        n_open_nodes == 0 ||
        n_cover_edges == 0 ||
        cover_edge_nodes == NULL ||
        n_node_map == 0 ||
        node_map == NULL)
    {
        return -1;
    }


    aperture_coupling_free(
        coupling
    );


    coupling->n_open_nodes =
    n_open_nodes;

    coupling->n_cover_edges =
    n_cover_edges;

    coupling->n_node_map =
    n_node_map;


    coupling->node_map = calloc(
        n_node_map,
        sizeof(ApertureNodeMapEntry)
    );

    coupling->edge_nodes = calloc(
        2 * n_cover_edges,
        sizeof(size_t)
    );

    coupling->edge_class = calloc(
        n_cover_edges,
        sizeof(ApertureEdgeClass)
    );

    coupling->B_gamma = calloc(
        n_open_nodes * n_cover_edges,
        sizeof(double)
    );


    if (coupling->node_map == NULL ||
        coupling->edge_nodes == NULL ||
        coupling->edge_class == NULL ||
        coupling->B_gamma == NULL)
    {
        aperture_coupling_free(
            coupling
        );

        return -1;
    }


    memcpy(
        coupling->node_map,
        node_map,
        n_node_map
        * sizeof(ApertureNodeMapEntry)
    );


    memcpy(
        coupling->edge_nodes,
        cover_edge_nodes,
        2
        * n_cover_edges
        * sizeof(size_t)
    );


    /*
     * Проверяем open node ids.
     */
    for (size_t i = 0; i < n_node_map; ++i)
    {
        if (coupling->node_map[i].open_node >=
            n_open_nodes)
        {
            fprintf(
                stderr,
                "ERROR: aperture open node %zu outside open mesh.\n",
                coupling->node_map[i].open_node
            );

            aperture_coupling_free(
                coupling
            );

            return -1;
        }
    }


    if (aperture_coupling_classify_edges(
        coupling
    ) != 0)
    {
        aperture_coupling_free(
            coupling
        );

        return -1;
    }


    if (aperture_coupling_build_B_gamma(
        coupling
    ) != 0)
    {
        aperture_coupling_free(
            coupling
        );

        return -1;
    }


    coupling->initialized =
    1;


    return 0;
}


/*
 * ============================================================
 * I_GAMMA = B_GAMMA I_COVER
 * ============================================================
 */
int aperture_coupling_compute_node_source(
    const ApertureCoupling *coupling,
    const double *cover_current,
    double *node_source
)
{
    if (coupling == NULL ||
        !coupling->initialized ||
        coupling->B_gamma == NULL ||
        cover_current == NULL ||
        node_source == NULL)
    {
        return -1;
    }


    size_t nv =
    coupling->n_open_nodes;

    size_t ne =
    coupling->n_cover_edges;


    for (size_t j = 0; j < nv; ++j)
    {
        double value =
        0.0;


        for (size_t e = 0; e < ne; ++e)
        {
            value +=
            coupling->B_gamma[
                matrix_index(
                    j,
                    e,
                    ne
                )
            ]
            *
            cover_current[e];
        }


        node_source[j] =
        value;
    }


    return 0;
}


/*
 * ============================================================
 * NET CURRENT
 * ============================================================
 */
double aperture_coupling_compute_net_current(
    const ApertureCoupling *coupling,
    const double *cover_current
)
{
    if (coupling == NULL ||
        !coupling->initialized ||
        coupling->B_gamma == NULL ||
        cover_current == NULL)
    {
        return NAN;
    }


    size_t nv =
    coupling->n_open_nodes;

    size_t ne =
    coupling->n_cover_edges;


    double net_current =
    0.0;


    /*
     * Вместо явного построения I_gamma можно сразу:
     *
     *     sum_j B_gamma[j,e].
     */
    for (size_t e = 0; e < ne; ++e)
    {
        double column_sum =
        0.0;


        for (size_t j = 0; j < nv; ++j)
        {
            column_sum +=
            coupling->B_gamma[
                matrix_index(
                    j,
                    e,
                    ne
                )
            ];
        }


        net_current +=
        column_sum
        *
        cover_current[e];
    }


    return net_current;
}


/*
 * ============================================================
 * P_B_GAMMA = P B_GAMMA
 * ============================================================
 */
int aperture_coupling_build_P_B_gamma(
    ApertureCoupling *coupling,
    const double *P_open
)
{
    if (coupling == NULL ||
        !coupling->initialized ||
        coupling->B_gamma == NULL ||
        P_open == NULL)
    {
        return -1;
    }


    size_t nv =
    coupling->n_open_nodes;

    size_t ne =
    coupling->n_cover_edges;


    free(
        coupling->P_B_gamma
    );

    coupling->P_B_gamma = calloc(
        nv * ne,
        sizeof(double)
    );


    if (coupling->P_B_gamma == NULL)
    {
        coupling->has_P_B_gamma =
        0;

        return -1;
    }


    /*
     * --------------------------------------------------------
     * P_B_gamma[j,e]
     *
     *     =
     *
     * sum_m
     *
     *     P[j,m] B_gamma[m,e].
     * --------------------------------------------------------
     */
    for (size_t j = 0; j < nv; ++j)
    {
        for (size_t e = 0; e < ne; ++e)
        {
            double value =
            0.0;


            for (size_t m = 0; m < nv; ++m)
            {
                value +=
                P_open[
                    matrix_index(
                        j,
                        m,
                        nv
                    )
                ]
                *
                coupling->B_gamma[
                    matrix_index(
                        m,
                        e,
                        ne
                    )
                ];
            }


            coupling->P_B_gamma[
                matrix_index(
                    j,
                    e,
                    ne
                )
            ] =
            value;
        }
    }


    coupling->has_P_B_gamma =
    1;


    return 0;
}


/*
 * ============================================================
 * P I_GAMMA
 *
 *     =
 *
 * P_B_GAMMA I_COVER
 * ============================================================
 */
int aperture_coupling_compute_P_node_source(
    const ApertureCoupling *coupling,
    const double *cover_current,
    double *P_node_source
)
{
    if (coupling == NULL ||
        !coupling->initialized ||
        !coupling->has_P_B_gamma ||
        coupling->P_B_gamma == NULL ||
        cover_current == NULL ||
        P_node_source == NULL)
    {
        return -1;
    }


    size_t nv =
    coupling->n_open_nodes;

    size_t ne =
    coupling->n_cover_edges;


    for (size_t j = 0; j < nv; ++j)
    {
        double value =
        0.0;


        for (size_t e = 0; e < ne; ++e)
        {
            value +=
            coupling->P_B_gamma[
                matrix_index(
                    j,
                    e,
                    ne
                )
            ]
            *
            cover_current[e];
        }


        P_node_source[j] =
        value;
    }


    return 0;
}


/*
 * ============================================================
 * PRINT INFO
 * ============================================================
 */
void aperture_coupling_print_info(
    const ApertureCoupling *coupling
)
{
    if (coupling == NULL)
    {
        return;
    }


    printf("\n");
    printf("============================================================\n");
    printf("APERTURE COUPLING\n");
    printf("============================================================\n");


    printf(
        "Open mesh nodes     : %zu\n",
        coupling->n_open_nodes
    );


    printf(
        "PEC cover edges     : %zu\n",
        coupling->n_cover_edges
    );


    printf(
        "Boundary node map   : %zu\n",
        coupling->n_node_map
    );


    printf("\n");


    printf(
        "Perimeter edges     : %zu\n",
        coupling->n_perimeter_edges
    );


    printf(
        "Cut edges           : %zu\n",
        coupling->n_cut_edges
    );


    printf(
        "Internal edges      : %zu\n",
        coupling->n_internal_edges
    );


    printf("\n");


    printf(
        "B_gamma shape       : %zu x %zu\n",
        coupling->n_open_nodes,
        coupling->n_cover_edges
    );


    if (coupling->has_P_B_gamma)
    {
        printf(
            "P_B_gamma           : READY\n"
        );
    }
    else
    {
        printf(
            "P_B_gamma           : NOT BUILT\n"
        );
    }


    printf("\n");


    printf(
        "Transfer equation   : I_gamma = B_gamma I_cover\n"
    );


    printf("============================================================\n");
}
