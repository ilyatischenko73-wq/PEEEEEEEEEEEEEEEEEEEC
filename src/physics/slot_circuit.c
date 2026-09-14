#include "slot_circuit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


static int positive_finite(
    double value
)
{
    return isfinite(value)
    &&
    value > 0.0;
}


void slot_circuit_init(
    SlotCircuit *circuit
)
{
    if (circuit == NULL)
    {
        return;
    }

    memset(
        circuit,
        0,
        sizeof(*circuit)
    );
}


void slot_circuit_free(
    SlotCircuit *circuit
)
{
    if (circuit == NULL)
    {
        return;
    }

    free(
        circuit->cells
    );

    memset(
        circuit,
        0,
        sizeof(*circuit)
    );
}


int slot_line_model_build(
    SlotLineModel *model,
    double slot_width,
    double wall_span,
    double wall_thickness,
    double epsilon,
    double mu
)
{
    if (model == NULL ||
        !positive_finite(slot_width) ||
        !positive_finite(wall_span) ||
        !isfinite(wall_thickness) ||
        wall_thickness < 0.0 ||
        !positive_finite(epsilon) ||
        !positive_finite(mu))
    {
        return -1;
    }

    memset(
        model,
        0,
        sizeof(*model)
    );

    /*
     * Effective-width correction used in the classical
     * enclosure-aperture equivalent-circuit model.
     *
     * For an ideal zero-thickness PEC sheet:
     *
     *     w_e = w.
     */
    double effective_width =
    slot_width;

    if (wall_thickness > 0.0)
    {
        double argument =
        4.0
        * M_PI
        * slot_width
        / wall_thickness;

        if (!positive_finite(
            argument
        ))
        {
            return -1;
        }

        effective_width =
        slot_width
        -
        (
            5.0
            * wall_thickness
            /
            (
                4.0
                * M_PI
            )
        )
        *
        (
            1.0
            +
            log(
                argument
            )
        );
    }

    if (!positive_finite(
        effective_width
    ))
    {
        fprintf(
            stderr,
            "ERROR: slot effective width <= 0.\n"
        );

        return -1;
    }

    /*
     * Narrow-aperture domain for the closed-form approximation.
     */
    if (effective_width >=
        wall_span / sqrt(2.0))
    {
        fprintf(
            stderr,
            "ERROR: slot aperture is outside the narrow-slot "
            "range: we/b = %.9e.\n",
            effective_width / wall_span
        );

        return -1;
    }

    double ratio =
    effective_width
    / wall_span;

    double radicand =
    1.0
    -
    ratio
    * ratio;

    if (radicand <= 0.0)
    {
        return -1;
    }

    /*
     * q = fourth_root(1 - (we/b)^2)
     */
    double q =
    sqrt(
        sqrt(
            radicand
        )
    );

    double logarithm_argument =
    2.0
    *
    (
        1.0 + q
    )
    /
    (
        1.0 - q
    );

    if (!positive_finite(
        logarithm_argument
    ))
    {
        return -1;
    }

    double denominator =
    log(
        logarithm_argument
    );

    if (!positive_finite(
        denominator
    ))
    {
        return -1;
    }

    /*
     * Free-space Robinson/Gupta approximation:
     *
     *     Z0s =
     *     120*pi^2 /
     *     ln[
     *         2(1+q)/(1-q)
     *       ].
     *
     * For a homogeneous material medium we scale by eta/eta0.
     */
    const double mu0 =
    4.0e-7
    * M_PI;

    const double epsilon0 =
    8.8541878128e-12;

    double eta0 =
    sqrt(
        mu0
        / epsilon0
    );

    double eta =
    sqrt(
        mu
        / epsilon
    );

    double Z0 =
    (
        120.0
        * M_PI
        * M_PI
        / denominator
    )
    *
    eta
    / eta0;

    double velocity =
    1.0
    /
    sqrt(
        mu
        * epsilon
    );

    double L_per_length =
    Z0
    / velocity;

    double C_per_length =
    1.0
    /
    (
        Z0
        * velocity
    );

    if (!positive_finite(Z0) ||
        !positive_finite(velocity) ||
        !positive_finite(L_per_length) ||
        !positive_finite(C_per_length))
    {
        return -1;
    }

    model->slot_width =
    slot_width;

    model->wall_span =
    wall_span;

    model->wall_thickness =
    wall_thickness;

    model->epsilon =
    epsilon;

    model->mu =
    mu;

    model->effective_width =
    effective_width;

    model->characteristic_impedance =
    Z0;

    model->wave_velocity =
    velocity;

    model->inductance_per_length =
    L_per_length;

    model->capacitance_per_length =
    C_per_length;

    return 0;
}


int slot_circuit_load_cells(
    SlotCircuit *circuit,
    const char *filename,
    const SlotLineModel *model,
    size_t n_mesh_nodes
)
{
    if (circuit == NULL ||
        filename == NULL ||
        model == NULL ||
        n_mesh_nodes == 0)
    {
        return -1;
    }

    FILE *file =
    fopen(
        filename,
        "r"
    );

    if (file == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot open slot-cell file '%s'.\n",
            filename
        );

        return -1;
    }

    slot_circuit_free(
        circuit
    );

    circuit->line =
    *model;

    size_t capacity =
    16;

    circuit->cells =
    malloc(
        capacity
        * sizeof(SlotCell)
    );

    if (circuit->cells == NULL)
    {
        fclose(
            file
        );

        return -1;
    }

    char line[4096];

    while (fgets(
        line,
        sizeof(line),
        file
    ) != NULL)
    {
        char *p =
        line;

        while (*p == ' ' ||
               *p == '\t')
        {
            ++p;
        }

        if (*p == '\0' ||
            *p == '\n' ||
            *p == '\r' ||
            *p == '#')
        {
            continue;
        }

        size_t node_a =
        0;

        size_t node_b =
        0;

        double support_length =
        0.0;

        if (sscanf(
            p,
            "%zu %zu %lf",
            &node_a,
            &node_b,
            &support_length
        ) != 3)
        {
            fprintf(
                stderr,
                "ERROR: invalid slot-cell line: %s",
                line
            );

            fclose(
                file
            );

            slot_circuit_free(
                circuit
            );

            return -1;
        }

        if (node_a >= n_mesh_nodes ||
            node_b >= n_mesh_nodes ||
            node_a == node_b ||
            !positive_finite(support_length))
        {
            fprintf(
                stderr,
                "ERROR: invalid slot cell A=%zu B=%zu dl=%.9e.\n",
                node_a,
                node_b,
                support_length
            );

            fclose(
                file
            );

            slot_circuit_free(
                circuit
            );

            return -1;
        }

        if (circuit->n_cells ==
            capacity)
        {
            capacity *=
            2;

            SlotCell *new_cells =
            realloc(
                circuit->cells,
                capacity
                * sizeof(SlotCell)
            );

            if (new_cells == NULL)
            {
                fclose(
                    file
                );

                slot_circuit_free(
                    circuit
                );

                return -1;
            }

            circuit->cells =
            new_cells;
        }

        SlotCell *cell =
        &circuit->cells[
            circuit->n_cells
        ];

        cell->node_a =
        node_a;

        cell->node_b =
        node_b;

        cell->support_length =
        support_length;

        cell->inductance =
        model->inductance_per_length
        * support_length;

        cell->capacitance =
        model->capacitance_per_length
        * support_length;

        ++circuit->n_cells;
    }

    fclose(
        file
    );

    if (circuit->n_cells == 0)
    {
        fprintf(
            stderr,
            "ERROR: slot-cell file contains no cells.\n"
        );

        slot_circuit_free(
            circuit
        );

        return -1;
    }

    return 0;
}


void slot_circuit_print_info(
    const SlotCircuit *circuit
)
{
    if (circuit == NULL)
    {
        return;
    }

    printf("\n");
    printf("============================================================\n");
    printf("DISTRIBUTED EQUIVALENT SLOT CELLS\n");
    printf("============================================================\n");

    printf(
        "Cells              : %zu\n",
        circuit->n_cells
    );

    printf(
        "Slot width         : %.9e m\n",
        circuit->line.slot_width
    );

    printf(
        "Effective width    : %.9e m\n",
        circuit->line.effective_width
    );

    printf(
        "Wall span          : %.9e m\n",
        circuit->line.wall_span
    );

    printf(
        "Wall thickness     : %.9e m\n",
        circuit->line.wall_thickness
    );

    printf(
        "Z0 slot            : %.9e Ohm\n",
        circuit->line.characteristic_impedance
    );

    printf(
        "Wave velocity      : %.9e m/s\n",
        circuit->line.wave_velocity
    );

    printf(
        "L'                 : %.9e H/m\n",
        circuit->line.inductance_per_length
    );

    printf(
        "C'                 : %.9e F/m\n",
        circuit->line.capacitance_per_length
    );

    printf("\n");
    printf(
        " cell      A      B        dl [m]            L [H]            C [F]\n"
    );
    printf(
        "------------------------------------------------------------------------\n"
    );

    for (size_t m = 0;
         m < circuit->n_cells;
    ++m)
    {
        const SlotCell *cell =
        &circuit->cells[m];

        printf(
            "%5zu  %5zu  %5zu  %.9e  %.9e  %.9e\n",
            m,
            cell->node_a,
            cell->node_b,
            cell->support_length,
            cell->inductance,
            cell->capacitance
        );
    }

    printf("============================================================\n");
}
