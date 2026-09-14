#include "lightning.h"

#include <math.h>
#include <stdio.h>
#include <string.h>


/*
 * ============================================================
 * DEFAULT
 * ============================================================
 */
void lightning_source_default(
    LightningSource *source
)
{
    if (source == NULL)
    {
        return;
    }

    source->wave_type =
    LIGHTNING_WAVE_DOUBLE_EXPONENTIAL;

    /*
     * Стартовая модель тока:
     *
     *     I(t) =
     *
     *     102900 [
     *         exp(-1500 t)
     *         -
     *         exp(-1e6 t)
     *     ].
     */
    source->K =
    102900.0;

    source->alpha =
    1500.0;

    source->beta =
    1.0e6;

    source->delay =
    0.0;

    /*
     * Это только безопасные начальные значения.
     *
     * Перед реальным расчётом их нужно заменить.
     */
    source->strike_node =
    0;

    source->return_node =
    0;
}


/*
 * ============================================================
 * VALIDATE
 * ============================================================
 */
int lightning_source_validate(
    const LightningSource *source,
    size_t n_nodes
)
{
    if (source == NULL ||
        n_nodes == 0)
    {
        return -1;
    }

    if (source->wave_type !=
        LIGHTNING_WAVE_DOUBLE_EXPONENTIAL)
    {
        fprintf(
            stderr,
            "ERROR: unsupported lightning waveform.\n"
        );

        return -1;
    }

    if (!isfinite(source->K) ||
        source->K <= 0.0)
    {
        fprintf(
            stderr,
            "ERROR: lightning K must be positive.\n"
        );

        return -1;
    }

    if (!isfinite(source->alpha) ||
        !isfinite(source->beta) ||
        source->alpha <= 0.0 ||
        source->beta <= 0.0)
    {
        fprintf(
            stderr,
            "ERROR: lightning alpha and beta must be positive.\n"
        );

        return -1;
    }

    /*
     * Для положительного импульса:
     *
     *     exp(-alpha*t) - exp(-beta*t) > 0
     *
     * при:
     *
     *     beta > alpha.
     */
    if (source->beta <= source->alpha)
    {
        fprintf(
            stderr,
            "ERROR: lightning requires beta > alpha.\n"
        );

        return -1;
    }

    if (!isfinite(source->delay) ||
        source->delay < 0.0)
    {
        fprintf(
            stderr,
            "ERROR: lightning delay must be non-negative.\n"
        );

        return -1;
    }

    if (source->strike_node >= n_nodes)
    {
        fprintf(
            stderr,
            "ERROR: strike node %zu is outside mesh.\n",
            source->strike_node
        );

        return -1;
    }

    if (source->return_node >= n_nodes)
    {
        fprintf(
            stderr,
            "ERROR: return node %zu is outside mesh.\n",
            source->return_node
        );

        return -1;
    }

    if (source->strike_node ==
        source->return_node)
    {
        fprintf(
            stderr,
            "ERROR: strike and return nodes must be different.\n"
        );

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * I_L(t)
 * ============================================================
 */
double lightning_current(
    const LightningSource *source,
    double time
)
{
    if (source == NULL ||
        !isfinite(time))
    {
        return NAN;
    }

    if (time < source->delay)
    {
        return 0.0;
    }

    double tau =
    time - source->delay;

    if (source->wave_type ==
        LIGHTNING_WAVE_DOUBLE_EXPONENTIAL)
    {
        return source->K
        * (
            exp(-source->alpha * tau)
            -
            exp(-source->beta * tau)
        );
    }

    return NAN;
}


/*
 * ============================================================
 * ПИК ДВОЙНОЙ ЭКСПОНЕНТЫ
 * ============================================================
 *
 *     I(tau) =
 *
 *     K [
 *         exp(-alpha*tau)
 *         -
 *         exp(-beta*tau)
 *       ].
 *
 *
 * Производная:
 *
 *     dI/dtau =
 *
 *     K [
 *         -alpha exp(-alpha*tau)
 *         +
 *          beta exp(-beta*tau)
 *       ].
 *
 *
 * При dI/dtau = 0:
 *
 *     beta exp(-beta*tau)
 *
 *     =
 *
 *     alpha exp(-alpha*tau).
 *
 *
 * Отсюда:
 *
 *     tau_peak =
 *
 *     ln(beta/alpha)
 *     ----------------
 *       beta-alpha.
 */
int lightning_peak(
    const LightningSource *source,
    double *time_peak,
    double *current_peak
)
{
    if (source == NULL ||
        time_peak == NULL ||
        current_peak == NULL)
    {
        return -1;
    }

    if (source->beta <= source->alpha ||
        source->alpha <= 0.0)
    {
        return -1;
    }

    double tau_peak =
    log(
        source->beta
        / source->alpha
    )
    / (
        source->beta
        - source->alpha
    );

    *time_peak =
    source->delay
    + tau_peak;

    *current_peak =
    source->K
    * (
        exp(-source->alpha * tau_peak)
        -
        exp(-source->beta * tau_peak)
    );

    return 0;
}


/*
 * ============================================================
 * NODE CURRENT VECTOR
 * ============================================================
 */
int lightning_build_node_current(
    const LightningSource *source,
    size_t n_nodes,
    double time,
    double *node_current
)
{
    if (source == NULL ||
        node_current == NULL ||
        !isfinite(time))
    {
        return -1;
    }

    if (lightning_source_validate(
        source,
        n_nodes
    ) != 0)
    {
        return -1;
    }

    /*
     * Сначала:
     *
     *     I_s = 0.
     */
    memset(
        node_current,
        0,
        n_nodes * sizeof(double)
    );

    double current =
    lightning_current(
        source,
        time
    );

    if (!isfinite(current))
    {
        return -1;
    }

    /*
     * Ввод тока в корпус.
     */
    node_current[source->strike_node] =
    current;

    /*
     * Возврат тока.
     */
    node_current[source->return_node] =
    -current;

    return 0;
}


/*
 * ============================================================
 * NEAREST NODE
 * ============================================================
 */
int lightning_find_nearest_node(
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
 * PRINT
 * ============================================================
 */
void lightning_print_info(
    const LightningSource *source,
    const Mesh *mesh
)
{
    if (source == NULL)
    {
        return;
    }

    double time_peak =
    NAN;

    double current_peak =
    NAN;

    lightning_peak(
        source,
        &time_peak,
        &current_peak
    );


    printf("\n");
    printf("============================================================\n");
    printf("LIGHTNING CURRENT SOURCE\n");
    printf("============================================================\n");

    printf(
        "Waveform          : DOUBLE EXPONENTIAL\n"
    );

    printf(
        "K                 : %.9e A\n",
        source->K
    );

    printf(
        "alpha             : %.9e 1/s\n",
        source->alpha
    );

    printf(
        "beta              : %.9e 1/s\n",
        source->beta
    );

    printf(
        "delay             : %.9e s\n",
        source->delay
    );

    printf(
        "Peak time         : %.9e s\n",
        time_peak
    );

    printf(
        "Peak current      : %.9e A\n",
        current_peak
    );


    printf("\n");

    printf(
        "Strike node       : %zu\n",
        source->strike_node
    );

    printf(
        "Return node       : %zu\n",
        source->return_node
    );


    /*
     * Если сетка доступна, дополнительно печатаем
     * координаты выбранных узлов.
     */
    if (mesh != NULL &&
        mesh->xyz != NULL &&
        source->strike_node < mesh->n_nodes &&
        source->return_node < mesh->n_nodes)
    {
        const double *strike =
        &mesh->xyz[
            3 * source->strike_node
        ];

        const double *return_point =
        &mesh->xyz[
            3 * source->return_node
        ];


        printf(
            "Strike position   : %.9e %.9e %.9e m\n",
            strike[0],
            strike[1],
            strike[2]
        );

        printf(
            "Return position   : %.9e %.9e %.9e m\n",
            return_point[0],
            return_point[1],
            return_point[2]
        );
    }


    printf("============================================================\n");
}
