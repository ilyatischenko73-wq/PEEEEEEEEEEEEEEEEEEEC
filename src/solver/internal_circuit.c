#include "internal_circuit.h"

#include <math.h>
#include <stdio.h>
#include <string.h>


/*
 * ============================================================
 * DEFAULT
 * ============================================================
 */
void internal_circuit_shunt_default(
    InternalCircuitShunt *shunt
)
{
    if (shunt == NULL)
    {
        return;
    }

    memset(
        shunt,
        0,
        sizeof(*shunt)
    );

    shunt->type =
    INTERNAL_CIRCUIT_RESISTIVE_SHUNT;

    /*
     * То же стартовое значение, что использовалось
     * в Python Stage 2.
     */
    shunt->resistance =
    1000.0;

    shunt->node_a =
    0;

    shunt->node_b =
    0;

    shunt->mapping_error_a =
    NAN;

    shunt->mapping_error_b =
    NAN;

    shunt->nodes_mapped =
    0;

    shunt->initialized =
    1;
}


/*
 * ============================================================
 * ПОИСК БЛИЖАЙШЕГО УЗЛА
 * ============================================================
 */
int internal_circuit_find_nearest_node(
    const Mesh *mesh,
    const double point[3],
    size_t *node_id,
    double *distance
)
{
    if (mesh == NULL ||
        point == NULL ||
        node_id == NULL ||
        mesh->xyz == NULL ||
        mesh->n_nodes == 0)
    {
        return -1;
    }


    size_t best_node =
    0;


    double dx =
    mesh->xyz[0] - point[0];

    double dy =
    mesh->xyz[1] - point[1];

    double dz =
    mesh->xyz[2] - point[2];


    double best_distance_squared =
    dx * dx
    + dy * dy
    + dz * dz;


    for (size_t j = 1; j < mesh->n_nodes; ++j)
    {
        const double *r =
        &mesh->xyz[3 * j];


        dx =
        r[0] - point[0];

        dy =
        r[1] - point[1];

        dz =
        r[2] - point[2];


        double distance_squared =
        dx * dx
        + dy * dy
        + dz * dz;


        if (distance_squared <
            best_distance_squared)
        {
            best_distance_squared =
            distance_squared;

            best_node =
            j;
        }
    }


    *node_id =
    best_node;


    if (distance != NULL)
    {
        *distance =
        sqrt(
            best_distance_squared
        );
    }


    return 0;
}


/*
 * ============================================================
 * ЗАДАНИЕ УЗЛОВ НАПРЯМУЮ
 * ============================================================
 */
int internal_circuit_shunt_set_nodes(
    InternalCircuitShunt *shunt,
    const Mesh *mesh,
    size_t node_a,
    size_t node_b
)
{
    if (shunt == NULL ||
        mesh == NULL ||
        mesh->xyz == NULL ||
        node_a >= mesh->n_nodes ||
        node_b >= mesh->n_nodes)
    {
        return -1;
    }


    if (node_a == node_b)
    {
        fprintf(
            stderr,
            "ERROR: internal shunt terminals map to the same node.\n"
        );

        return -1;
    }


    shunt->node_a =
    node_a;

    shunt->node_b =
    node_b;


    const double *a =
    &mesh->xyz[3 * node_a];

    const double *b =
    &mesh->xyz[3 * node_b];


    for (int c = 0; c < 3; ++c)
    {
        shunt->requested_a[c] =
        a[c];

        shunt->requested_b[c] =
        b[c];

        shunt->mapped_a[c] =
        a[c];

        shunt->mapped_b[c] =
        b[c];
    }


    shunt->mapping_error_a =
    0.0;

    shunt->mapping_error_b =
    0.0;


    shunt->nodes_mapped =
    1;


    return 0;
}


/*
 * ============================================================
 * ЗАДАНИЕ ШУНТА ПО КООРДИНАТАМ
 * ============================================================
 */
int internal_circuit_shunt_map_coordinates(
    InternalCircuitShunt *shunt,
    const Mesh *mesh,
    const double point_a[3],
    const double point_b[3]
)
{
    if (shunt == NULL ||
        mesh == NULL ||
        point_a == NULL ||
        point_b == NULL)
    {
        return -1;
    }


    size_t node_a =
    0;

    size_t node_b =
    0;


    double error_a =
    0.0;

    double error_b =
    0.0;


    if (internal_circuit_find_nearest_node(
        mesh,
        point_a,
        &node_a,
        &error_a
    ) != 0)
    {
        return -1;
    }


    if (internal_circuit_find_nearest_node(
        mesh,
        point_b,
        &node_b,
        &error_b
    ) != 0)
    {
        return -1;
    }


    if (node_a == node_b)
    {
        fprintf(
            stderr,
            "ERROR: the two internal shunt terminals mapped "
            "to the same mesh node %zu.\n",
            node_a
        );

        return -1;
    }


    shunt->node_a =
    node_a;

    shunt->node_b =
    node_b;


    for (int c = 0; c < 3; ++c)
    {
        shunt->requested_a[c] =
        point_a[c];

        shunt->requested_b[c] =
        point_b[c];


        shunt->mapped_a[c] =
        mesh->xyz[
            3 * node_a + c
        ];

        shunt->mapped_b[c] =
        mesh->xyz[
            3 * node_b + c
        ];
    }


    shunt->mapping_error_a =
    error_a;

    shunt->mapping_error_b =
    error_b;


    shunt->nodes_mapped =
    1;


    return 0;
}


/*
 * ============================================================
 * VALIDATE
 * ============================================================
 */
int internal_circuit_shunt_validate(
    const InternalCircuitShunt *shunt,
    size_t n_nodes
)
{
    if (shunt == NULL ||
        !shunt->initialized ||
        n_nodes == 0)
    {
        return -1;
    }


    if (shunt->type !=
        INTERNAL_CIRCUIT_RESISTIVE_SHUNT)
    {
        fprintf(
            stderr,
            "ERROR: unsupported internal circuit type.\n"
        );

        return -1;
    }


    if (!shunt->nodes_mapped)
    {
        fprintf(
            stderr,
            "ERROR: internal shunt nodes have not been mapped.\n"
        );

        return -1;
    }


    if (shunt->node_a >= n_nodes ||
        shunt->node_b >= n_nodes)
    {
        fprintf(
            stderr,
            "ERROR: internal shunt node is outside mesh.\n"
        );

        return -1;
    }


    if (shunt->node_a ==
        shunt->node_b)
    {
        fprintf(
            stderr,
            "ERROR: internal shunt terminals must be different.\n"
        );

        return -1;
    }


    /*
     * R = 0 разрешён:
     *
     * это идеальное короткое замыкание.
     */
    if (!isfinite(shunt->resistance) ||
        shunt->resistance < 0.0)
    {
        fprintf(
            stderr,
            "ERROR: internal shunt resistance must be >= 0 Ohm.\n"
        );

        return -1;
    }


    return 0;
}


/*
 * ============================================================
 * b_R
 * ============================================================
 *
 * Положительный ток:
 *
 *          A --------> B
 *
 *
 * Поэтому:
 *
 *     b_R[A] = -1
 *     b_R[B] = +1.
 */
int internal_circuit_build_incidence_vector(
    const InternalCircuitShunt *shunt,
    size_t n_nodes,
    double *b_R
)
{
    if (shunt == NULL ||
        b_R == NULL)
    {
        return -1;
    }


    if (internal_circuit_shunt_validate(
        shunt,
        n_nodes
    ) != 0)
    {
        return -1;
    }


    memset(
        b_R,
        0,
        n_nodes * sizeof(double)
    );


    b_R[shunt->node_a] =
    -1.0;

    b_R[shunt->node_b] =
    +1.0;


    return 0;
}


/*
 * ============================================================
 * V_SHUNT
 * ============================================================
 *
 * Определяем:
 *
 *     V_shunt =
 *
 *     V_A - V_B.
 *
 *
 * При нашей ориентации:
 *
 *     b_R^T V =
 *
 *     V_B - V_A
 *
 * поэтому:
 *
 *     V_shunt =
 *
 *     -b_R^T V.
 */
double internal_circuit_shunt_voltage(
    const InternalCircuitShunt *shunt,
    const double *node_voltage,
    size_t n_nodes
)
{
    if (shunt == NULL ||
        node_voltage == NULL)
    {
        return NAN;
    }


    if (internal_circuit_shunt_validate(
        shunt,
        n_nodes
    ) != 0)
    {
        return NAN;
    }


    return node_voltage[shunt->node_a]
    - node_voltage[shunt->node_b];
}


/*
 * ============================================================
 * I_SHUNT ИЗ V
 * ============================================================
 *
 * Эта функция НЕ используется для формирования основной
 * Stage-2 системы.
 *
 * Она нужна только для диагностической проверки.
 */
double internal_circuit_shunt_current_from_voltage(
    const InternalCircuitShunt *shunt,
    const double *node_voltage,
    size_t n_nodes
)
{
    if (shunt == NULL ||
        node_voltage == NULL)
    {
        return NAN;
    }


    if (internal_circuit_shunt_validate(
        shunt,
        n_nodes
    ) != 0)
    {
        return NAN;
    }


    /*
     * При идеальном КЗ:
     *
     *     R = 0
     *
     * нельзя вычислять ток как V/R.
     *
     * Ток определяется всей PEEC-системой.
     */
    if (shunt->resistance == 0.0)
    {
        return NAN;
    }


    double voltage =
    internal_circuit_shunt_voltage(
        shunt,
        node_voltage,
        n_nodes
    );


    return voltage
    / shunt->resistance;
}


/*
 * ============================================================
 * НЕВЯЗКА ЗАКОНА ШУНТА
 * ============================================================
 *
 * Из уравнения:
 *
 *     b_R^T V + R I_shunt = 0
 *
 * следует:
 *
 *     V_A - V_B - R I_shunt = 0.
 *
 *
 * Возвращаем именно:
 *
 *     V_A - V_B - R I_shunt.
 */
double internal_circuit_shunt_residual(
    const InternalCircuitShunt *shunt,
    const double *node_voltage,
    size_t n_nodes,
    double shunt_current
)
{
    if (shunt == NULL ||
        node_voltage == NULL ||
        !isfinite(shunt_current))
    {
        return NAN;
    }


    if (internal_circuit_shunt_validate(
        shunt,
        n_nodes
    ) != 0)
    {
        return NAN;
    }


    double voltage =
    internal_circuit_shunt_voltage(
        shunt,
        node_voltage,
        n_nodes
    );


    return voltage
    - shunt->resistance
    * shunt_current;
}


/*
 * ============================================================
 * PRINT
 * ============================================================
 */
void internal_circuit_shunt_print_info(
    const InternalCircuitShunt *shunt
)
{
    if (shunt == NULL)
    {
        return;
    }


    printf("\n");
    printf("============================================================\n");
    printf("INTERNAL CIRCUIT / SHUNT\n");
    printf("============================================================\n");


    printf(
        "Type               : RESISTIVE SHUNT\n"
    );


    printf(
        "Positive current   : A -> B\n"
    );


    printf(
        "A node             : %zu\n",
        shunt->node_a
    );


    printf(
        "A requested xyz    : %.9e %.9e %.9e m\n",
        shunt->requested_a[0],
        shunt->requested_a[1],
        shunt->requested_a[2]
    );


    printf(
        "A mesh xyz         : %.9e %.9e %.9e m\n",
        shunt->mapped_a[0],
        shunt->mapped_a[1],
        shunt->mapped_a[2]
    );


    printf(
        "A mapping error    : %.9e m\n",
        shunt->mapping_error_a
    );


    printf("\n");


    printf(
        "B node             : %zu\n",
        shunt->node_b
    );


    printf(
        "B requested xyz    : %.9e %.9e %.9e m\n",
        shunt->requested_b[0],
        shunt->requested_b[1],
        shunt->requested_b[2]
    );


    printf(
        "B mesh xyz         : %.9e %.9e %.9e m\n",
        shunt->mapped_b[0],
        shunt->mapped_b[1],
        shunt->mapped_b[2]
    );


    printf(
        "B mapping error    : %.9e m\n",
        shunt->mapping_error_b
    );


    printf("\n");


    printf(
        "Resistance         : %.9e Ohm\n",
        shunt->resistance
    );


    if (shunt->resistance == 0.0)
    {
        printf(
            "Shunt model        : IDEAL SHORT CIRCUIT\n"
        );
    }
    else
    {
        printf(
            "Shunt model        : RESISTIVE\n"
        );
    }


    printf("\n");


    printf(
        "Branch equation    : b_R^T V + R I_shunt = 0\n"
    );


    printf(
        "Voltage definition : V_shunt = V_A - V_B\n"
    );


    printf("============================================================\n");
}
