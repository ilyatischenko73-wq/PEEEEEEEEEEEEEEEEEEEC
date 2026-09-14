#include "transient.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static size_t transient_matrix_index(
    size_t row,
    size_t col,
    size_t n_cols
)
{
    return row * n_cols + col;
}


static double transient_scheme_theta(
    TransientScheme scheme
)
{
    if (scheme == TRANSIENT_BACKWARD_EULER)
    {
        return 1.0;
    }

    if (scheme == TRANSIENT_TRAPEZOIDAL)
    {
        return 0.5;
    }

    return NAN;
}


static int transient_scheme_is_valid(
    TransientScheme scheme
)
{
    return
    scheme == TRANSIENT_BACKWARD_EULER
    ||
    scheme == TRANSIENT_TRAPEZOIDAL;
}


void transient_system_init(
    TransientSystem *system
)
{
    if (system == NULL)
    {
        return;
    }

    memset(
        system,
        0,
        sizeof(*system)
    );

    matrix_init(
        &system->AS
    );

    lu_init(
        &system->lu
    );
}


void transient_system_free(
    TransientSystem *system
)
{
    if (system == NULL)
    {
        return;
    }

    free(
        system->F
    );

    matrix_free(
        &system->AS
    );

    free(
        system->b_R
    );

    free(
        system->shunt_voltage_operator
    );

    free(
        system->slot_mass_matrix
    );

    free(
        system->slot_voltage_operators
    );

    free(
        system->system_matrix
    );

    free(
        system->rhs
    );

    free(
        system->solution
    );

    lu_free(
        &system->lu
    );

    memset(
        system,
        0,
        sizeof(*system)
    );

    matrix_init(
        &system->AS
    );

    lu_init(
        &system->lu
    );
}


size_t transient_current_offset(
    const TransientSystem *system
)
{
    (void)system;

    return 0;
}


size_t transient_voltage_offset(
    const TransientSystem *system
)
{
    if (system == NULL)
    {
        return 0;
    }

    return system->n_edges;
}


size_t transient_slot_current_offset(
    const TransientSystem *system
)
{
    if (system == NULL ||
        system->stage != TRANSIENT_SLOT_CELLS)
    {
        return (size_t)-1;
    }

    return
    system->n_edges
    +
    system->n_nodes;
}


size_t transient_shunt_index(
    const TransientSystem *system
)
{
    if (system == NULL)
    {
        return (size_t)-1;
    }

    if (system->stage == TRANSIENT_STAGE_2)
    {
        return
        system->n_edges
        +
        system->n_nodes;
    }

    if (system->stage == TRANSIENT_SLOT_CELLS &&
        system->shunt != NULL)
    {
        return
        transient_slot_current_offset(
            system
        )
        +
        system->n_slot_cells;
    }

    return (size_t)-1;
}


static int transient_validate_common_input(
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const PotentialMatrix *P,
    double dt
)
{
    if (L == NULL ||
        R == NULL ||
        A == NULL ||
        P == NULL)
    {
        fprintf(
            stderr,
            "ERROR: transient PEEC matrix pointer is NULL.\n"
        );

        return -1;
    }

    if (!isfinite(dt) ||
        dt <= 0.0)
    {
        fprintf(
            stderr,
            "ERROR: transient dt must be positive.\n"
        );

        return -1;
    }

    if (L->data == NULL ||
        R->diagonal == NULL ||
        A->data == NULL ||
        P->data == NULL)
    {
        fprintf(
            stderr,
            "ERROR: transient PEEC matrix data is NULL.\n"
        );

        return -1;
    }

    if (L->n == 0 ||
        P->n == 0)
    {
        return -1;
    }

    if (R->n != L->n)
    {
        fprintf(
            stderr,
            "ERROR: R dimension does not match L dimension.\n"
        );

        return -1;
    }

    if (A->rows != L->n ||
        A->cols != P->n)
    {
        fprintf(
            stderr,
            "ERROR: A dimensions are inconsistent with L and P.\n"
        );

        return -1;
    }

    return 0;
}


static int transient_allocate_common(
    TransientSystem *system
)
{
    if (system == NULL ||
        system->n_total == 0 ||
        system->n_nodes == 0)
    {
        return -1;
    }

    system->F = calloc(
        system->n_nodes,
        sizeof(double)
    );

    system->system_matrix = calloc(
        system->n_total * system->n_total,
        sizeof(double)
    );

    system->rhs = calloc(
        system->n_total,
        sizeof(double)
    );

    system->solution = calloc(
        system->n_total,
        sizeof(double)
    );

    if (system->F == NULL ||
        system->system_matrix == NULL ||
        system->rhs == NULL ||
        system->solution == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate transient system storage.\n"
        );

        return -1;
    }

    return 0;
}


static int transient_build_electrostatic_operators(
    TransientSystem *system
)
{
    if (system == NULL ||
        system->A == NULL ||
        system->P == NULL ||
        system->F == NULL)
    {
        return -1;
    }

    if (potential_build_self_capacitance(
        system->P,
        system->F
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot build transient F operator.\n"
        );

        return -1;
    }

    /*
     * AS строится один раз.
     *
     * Здесь используем один поток, чтобы API создания transient
     * системы не зависел от OpenMP-настроек. Основная стоимость
     * расчета все равно приходится на L/P и LU.
     */
    if (incidence_build_AS(
        system->A,
        system->P,
        &system->AS,
        1
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot build transient A S operator.\n"
        );

        return -1;
    }

    return 0;
}


int transient_system_create_stage1(
    TransientSystem *system,
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const PotentialMatrix *P,
    TransientScheme scheme,
    double dt
)
{
    if (system == NULL ||
        !transient_scheme_is_valid(
            scheme
        ))
    {
        return -1;
    }

    if (transient_validate_common_input(
        L,
        R,
        A,
        P,
        dt
    ) != 0)
    {
        return -1;
    }

    transient_system_free(
        system
    );

    system->stage =
    TRANSIENT_STAGE_1;

    system->scheme =
    scheme;

    system->n_edges =
    L->n;

    system->n_nodes =
    P->n;

    system->n_total =
    system->n_edges
    +
    system->n_nodes;

    system->dt =
    dt;

    system->L =
    L;

    system->R =
    R;

    system->A =
    A;

    system->P =
    P;

    system->shunt =
    NULL;

    if (transient_allocate_common(
        system
    ) != 0 ||
    transient_build_electrostatic_operators(
        system
    ) != 0 ||
    transient_build_system_matrix(
        system
    ) != 0)
    {
        transient_system_free(
            system
        );

        return -1;
    }

    system->initialized =
    1;

    system->factorized =
    0;

    return 0;
}


int transient_system_create_stage2(
    TransientSystem *system,
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const PotentialMatrix *P,
    const InternalCircuitShunt *shunt,
    TransientScheme scheme,
    double dt
)
{
    if (system == NULL ||
        shunt == NULL ||
        !transient_scheme_is_valid(
            scheme
        ))
    {
        return -1;
    }

    if (transient_validate_common_input(
        L,
        R,
        A,
        P,
        dt
    ) != 0)
    {
        return -1;
    }

    if (internal_circuit_shunt_validate(
        shunt,
        P->n
    ) != 0)
    {
        return -1;
    }

    transient_system_free(
        system
    );

    system->stage =
    TRANSIENT_STAGE_2;

    system->scheme =
    scheme;

    system->n_edges =
    L->n;

    system->n_nodes =
    P->n;

    system->n_total =
    system->n_edges
    +
    system->n_nodes
    +
    1;

    system->dt =
    dt;

    system->L =
    L;

    system->R =
    R;

    system->A =
    A;

    system->P =
    P;

    system->shunt =
    shunt;

    system->b_R = calloc(
        system->n_nodes,
        sizeof(double)
    );

    system->shunt_voltage_operator = calloc(
        system->n_nodes,
        sizeof(double)
    );

    if (system->b_R == NULL ||
        system->shunt_voltage_operator == NULL)
    {
        transient_system_free(
            system
        );

        return -1;
    }

    if (transient_allocate_common(
        system
    ) != 0 ||
    transient_build_electrostatic_operators(
        system
    ) != 0 ||
    transient_build_shunt_operators(
        system
    ) != 0 ||
    transient_build_system_matrix(
        system
    ) != 0)
    {
        transient_system_free(
            system
        );

        return -1;
    }

    system->initialized =
    1;

    system->factorized =
    0;

    return 0;
}



static int transient_build_slot_operators(
    TransientSystem *system
)
{
    if (system == NULL ||
        system->stage != TRANSIENT_SLOT_CELLS ||
        system->P == NULL ||
        system->P->data == NULL ||
        system->slot_circuit == NULL ||
        system->slot_circuit->cells == NULL ||
        system->n_slot_cells == 0 ||
        system->slot_mass_matrix == NULL ||
        system->slot_voltage_operators == NULL)
    {
        return -1;
    }

    size_t nv =
    system->n_nodes;

    memset(
        system->slot_mass_matrix,
        0,
        nv * nv * sizeof(double)
    );

    memset(
        system->slot_voltage_operators,
        0,
        system->n_slot_cells
        * nv
        * sizeof(double)
    );

    for (size_t m = 0;
         m < system->n_slot_cells;
    ++m)
    {
        const SlotCell *cell =
        &system->slot_circuit->cells[m];

        if (cell->node_a >= nv ||
            cell->node_b >= nv ||
            cell->node_a == cell->node_b ||
            !isfinite(cell->inductance) ||
            cell->inductance <= 0.0 ||
            !isfinite(cell->capacitance) ||
            cell->capacitance <= 0.0)
        {
            return -1;
        }

        double *c =
        &system->slot_voltage_operators[
            m * nv
        ];

        /*
         * c_m^T = b_m^T S,
         *
         * b_m[A] = -1,
         * b_m[B] = +1.
         *
         * Therefore:
         *
         * c_m[a] =
         *
         *     S[B,a] - S[A,a]
         *
         *   = (P[B,a] - P[A,a]) / P[a,a].
         */
        for (size_t a = 0;
             a < nv;
        ++a)
        {
            double P_aa =
            system->P->data[
                a * nv + a
            ];

            if (!isfinite(P_aa) ||
                P_aa <= 0.0)
            {
                return -1;
            }

            c[a] =
            (
                system->P->data[
                    cell->node_b * nv + a
                ]
                -
                system->P->data[
                    cell->node_a * nv + a
                ]
            )
            / P_aa;
        }

        /*
         * M_slot += C_m b_m c_m^T.
         *
         * Only rows A and B are nonzero because b_m is sparse:
         *
         *     b_m[A] = -1
         *     b_m[B] = +1.
         */
        for (size_t a = 0;
             a < nv;
        ++a)
        {
            double value =
            cell->capacitance
            * c[a];

            system->slot_mass_matrix[
                cell->node_a * nv + a
            ] -=
            value;

            system->slot_mass_matrix[
                cell->node_b * nv + a
            ] +=
            value;
        }
    }

    return 0;
}


int transient_system_create_slot_cells(
    TransientSystem *system,
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const PotentialMatrix *P,
    const SlotCircuit *slot_circuit,
    const InternalCircuitShunt *shunt,
    TransientScheme scheme,
    double dt
)
{
    if (system == NULL ||
        slot_circuit == NULL ||
        slot_circuit->cells == NULL ||
        slot_circuit->n_cells == 0 ||
        !transient_scheme_is_valid(
            scheme
        ))
    {
        return -1;
    }

    if (transient_validate_common_input(
        L,
        R,
        A,
        P,
        dt
    ) != 0)
    {
        return -1;
    }

    if (shunt != NULL &&
        internal_circuit_shunt_validate(
            shunt,
            P->n
        ) != 0)
    {
        return -1;
    }

    transient_system_free(
        system
    );

    system->stage =
    TRANSIENT_SLOT_CELLS;

    system->scheme =
    scheme;

    system->n_edges =
    L->n;

    system->n_nodes =
    P->n;

    system->n_slot_cells =
    slot_circuit->n_cells;

    system->n_total =
    system->n_edges
    +
    system->n_nodes
    +
    system->n_slot_cells
    +
    (
        shunt != NULL
        ? 1
        : 0
    );

    system->dt =
    dt;

    system->L =
    L;

    system->R =
    R;

    system->A =
    A;

    system->P =
    P;

    system->slot_circuit =
    slot_circuit;

    system->shunt =
    shunt;

    system->slot_mass_matrix = calloc(
        system->n_nodes
        * system->n_nodes,
        sizeof(double)
    );

    system->slot_voltage_operators = calloc(
        system->n_slot_cells
        * system->n_nodes,
        sizeof(double)
    );

    if (system->slot_mass_matrix == NULL ||
        system->slot_voltage_operators == NULL)
    {
        transient_system_free(
            system
        );

        return -1;
    }

    if (shunt != NULL)
    {
        system->b_R = calloc(
            system->n_nodes,
            sizeof(double)
        );

        system->shunt_voltage_operator = calloc(
            system->n_nodes,
            sizeof(double)
        );

        if (system->b_R == NULL ||
            system->shunt_voltage_operator == NULL)
        {
            transient_system_free(
                system
            );

            return -1;
        }
    }

    if (transient_allocate_common(
        system
    ) != 0 ||
    transient_build_electrostatic_operators(
        system
    ) != 0 ||
    transient_build_slot_operators(
        system
    ) != 0)
    {
        transient_system_free(
            system
        );

        return -1;
    }

    if (shunt != NULL)
    {
        /*
         * Reuse the same shunt operators as Stage 2. Temporarily
         * build them here because transient_build_shunt_operators()
         * accepts the slot-cell stage as well.
         */
        if (transient_build_shunt_operators(
            system
        ) != 0)
        {
            transient_system_free(
                system
            );

            return -1;
        }
    }

    if (transient_build_system_matrix(
        system
    ) != 0)
    {
        transient_system_free(
            system
        );

        return -1;
    }

    system->initialized =
    1;

    system->factorized =
    0;

    return 0;
}


int transient_build_shunt_operators(
    TransientSystem *system
)
{
    if (system == NULL ||
        (
            system->stage != TRANSIENT_STAGE_2 &&
            system->stage != TRANSIENT_SLOT_CELLS
        ) ||
        system->P == NULL ||
        system->P->data == NULL ||
        system->shunt == NULL ||
        system->b_R == NULL ||
        system->shunt_voltage_operator == NULL)
    {
        return -1;
    }

    size_t nv =
    system->n_nodes;

    if (internal_circuit_build_incidence_vector(
        system->shunt,
        nv,
        system->b_R
    ) != 0)
    {
        return -1;
    }

    size_t node_a =
    system->shunt->node_a;

    size_t node_b =
    system->shunt->node_b;

    /*
     * c^T = b_R^T S.
     *
     * При:
     *
     *     b_R[A] = -1,
     *     b_R[B] = +1,
     *
     * имеем:
     *
     *     c[a]
     *
     *     =
     *
     *     S[B,a] - S[A,a]
     *
     *     =
     *
     *     (P[B,a] - P[A,a]) / P[a,a].
     */
    for (size_t a = 0; a < nv; ++a)
    {
        double P_aa =
        system->P->data[a * nv + a];

        if (!isfinite(P_aa) ||
            P_aa <= 0.0)
        {
            return -1;
        }

        system->shunt_voltage_operator[a] =
        (
            system->P->data[node_b * nv + a]
            -
            system->P->data[node_a * nv + a]
        )
        / P_aa;
    }

    return 0;
}


int transient_build_system_matrix(
    TransientSystem *system
)
{
    if (system == NULL ||
        system->L == NULL ||
        system->R == NULL ||
        system->A == NULL ||
        system->P == NULL ||
        system->F == NULL ||
        system->AS.data == NULL ||
        system->system_matrix == NULL ||
        !transient_scheme_is_valid(
            system->scheme
        ) ||
        !isfinite(system->dt) ||
        system->dt <= 0.0)
    {
        return -1;
    }

    const double theta =
    transient_scheme_theta(
        system->scheme
    );

    size_t ne =
    system->n_edges;

    size_t nv =
    system->n_nodes;

    size_t nt =
    system->n_total;

    size_t i0 =
    transient_current_offset(
        system
    );

    size_t v0 =
    transient_voltage_offset(
        system
    );

    memset(
        system->system_matrix,
        0,
        nt * nt * sizeof(double)
    );


    /*
     * BLOCK (I,I):
     *
     *     L/dt + theta R.
     */
    for (size_t i = 0; i < ne; ++i)
    {
        for (size_t j = 0; j < ne; ++j)
        {
            double value =
            system->L->data[
                transient_matrix_index(
                    i,
                    j,
                    ne
                )
            ]
            / system->dt;

            if (i == j)
            {
                value +=
                theta
                * system->R->diagonal[i];
            }

            system->system_matrix[
                transient_matrix_index(
                    i0 + i,
                    i0 + j,
                    nt
                )
            ] =
            value;
        }
    }


    /*
     * BLOCK (I,V_c):
     *
     *     theta A S.
     */
    for (size_t e = 0; e < ne; ++e)
    {
        for (size_t j = 0; j < nv; ++j)
        {
            system->system_matrix[
                transient_matrix_index(
                    i0 + e,
                    v0 + j,
                    nt
                )
            ] =
            theta
            * system->AS.data[
                transient_matrix_index(
                    e,
                    j,
                    nv
                )
            ];
        }
    }


    /*
     * BLOCK (V_c,I):
     *
     *     -theta A^T.
     *
     *
     * A имеет размер Ne x Nv, поэтому:
     *
     *     A^T[j,e] = A[e,j].
     */
    for (size_t j = 0; j < nv; ++j)
    {
        for (size_t e = 0; e < ne; ++e)
        {
            system->system_matrix[
                transient_matrix_index(
                    v0 + j,
                    i0 + e,
                    nt
                )
            ] =
            -theta
            * system->A->data[
                transient_matrix_index(
                    e,
                    j,
                    nv
                )
            ];
        }
    }


    /*
     * BLOCK (V_c,V_c):
     *
     * Ordinary PEEC:
     *
     *     F/dt.
     *
     * Slot-cell model:
     *
     *     [diag(F) + M_slot] / dt,
     *
     * where:
     *
     *     M_slot = sum_m C_m b_m c_m^T.
     */
    for (size_t row = 0;
         row < nv;
    ++row)
    {
        for (size_t col = 0;
             col < nv;
        ++col)
        {
            double value =
            0.0;

            if (row == col)
            {
                value +=
                system->F[row];
            }

            if (system->stage ==
                TRANSIENT_SLOT_CELLS)
            {
                value +=
                system->slot_mass_matrix[
                    row * nv + col
                ];
            }

            system->system_matrix[
                transient_matrix_index(
                    v0 + row,
                    v0 + col,
                    nt
                )
            ] =
            value
            / system->dt;
        }
    }


    if (system->stage == TRANSIENT_SLOT_CELLS)
    {
        size_t s0 =
        transient_slot_current_offset(
            system
        );

        /*
         * NODE / SLOT INDUCTOR:
         *
         *     -theta B_slot I_L.
         *
         * Each sparse b_m has:
         *
         *     b_m[A] = -1
         *     b_m[B] = +1.
         */
        for (size_t m = 0;
             m < system->n_slot_cells;
        ++m)
        {
            const SlotCell *cell =
            &system->slot_circuit->cells[m];

            system->system_matrix[
                transient_matrix_index(
                    v0 + cell->node_a,
                    s0 + m,
                    nt
                )
            ] +=
            theta;

            system->system_matrix[
                transient_matrix_index(
                    v0 + cell->node_b,
                    s0 + m,
                    nt
                )
            ] -=
            theta;

            /*
             * SLOT INDUCTOR EQUATION:
             *
             *     L_m dI_m/dt
             *     + c_m^T V_c
             *     = 0.
             */
            const double *c =
            &system->slot_voltage_operators[
                m * nv
            ];

            for (size_t j = 0;
                 j < nv;
            ++j)
            {
                system->system_matrix[
                    transient_matrix_index(
                        s0 + m,
                        v0 + j,
                        nt
                    )
                ] =
                theta
                * c[j];
            }

            system->system_matrix[
                transient_matrix_index(
                    s0 + m,
                    s0 + m,
                    nt
                )
            ] =
            cell->inductance
            / system->dt;
        }

        if (system->shunt != NULL)
        {
            size_t s =
            transient_shunt_index(
                system
            );

            for (size_t j = 0;
                 j < nv;
            ++j)
            {
                system->system_matrix[
                    transient_matrix_index(
                        v0 + j,
                        s,
                        nt
                    )
                ] =
                -theta
                * system->b_R[j];

                system->system_matrix[
                    transient_matrix_index(
                        s,
                        v0 + j,
                        nt
                    )
                ] =
                system->shunt_voltage_operator[j];
            }

            system->system_matrix[
                transient_matrix_index(
                    s,
                    s,
                    nt
                )
            ] =
            system->shunt->resistance;
        }
    }


    if (system->stage == TRANSIENT_STAGE_2)
    {
        if (system->shunt == NULL ||
            system->b_R == NULL ||
            system->shunt_voltage_operator == NULL)
        {
            return -1;
        }

        size_t s =
        transient_shunt_index(
            system
        );

        /*
         * NODE / SHUNT:
         *
         *     -theta b_R I_shunt.
         */
        for (size_t j = 0; j < nv; ++j)
        {
            system->system_matrix[
                transient_matrix_index(
                    v0 + j,
                    s,
                    nt
                )
            ] =
            -theta
            * system->b_R[j];
        }


        /*
         * ALGEBRAIC SHUNT EQUATION:
         *
         *     b_R^T S V_c
         *     +
         *     R_shunt I_shunt
         *
         *     =
         *
         *     0.
         */
        for (size_t j = 0; j < nv; ++j)
        {
            system->system_matrix[
                transient_matrix_index(
                    s,
                    v0 + j,
                    nt
                )
            ] =
            system->shunt_voltage_operator[j];
        }

        system->system_matrix[
            transient_matrix_index(
                s,
                s,
                nt
            )
        ] =
        system->shunt->resistance;
    }

    system->factorized =
    0;

    return 0;
}


static int transient_build_rhs_common(
    TransientSystem *system,
    const double *edge_current_n,
    const double *node_self_voltage_n,
    const double *node_source_n,
    const double *node_source_np1,
    const double *edge_voltage_n,
    const double *edge_voltage_np1
)
{
    if (system == NULL ||
        !system->initialized ||
        edge_current_n == NULL ||
        node_self_voltage_n == NULL ||
        system->rhs == NULL ||
        system->F == NULL ||
        system->AS.data == NULL ||
        !transient_scheme_is_valid(
            system->scheme
        ))
    {
        return -1;
    }

    const double theta =
    transient_scheme_theta(
        system->scheme
    );

    const double one_minus_theta =
    1.0 - theta;

    size_t ne =
    system->n_edges;

    size_t nv =
    system->n_nodes;

    size_t nt =
    system->n_total;

    size_t i0 =
    transient_current_offset(
        system
    );

    size_t v0 =
    transient_voltage_offset(
        system
    );

    memset(
        system->rhs,
        0,
        nt * sizeof(double)
    );


    /*
     * EDGE EQUATION:
     *
     *     L dI/dt
     *     + R I
     *     + AS V_c
     *     + U_e
     *
     *     =
     *
     *     0.
     *
     *
     * Theta-схема:
     *
     *     (L/dt + theta R) I^(n+1)
     *     + theta AS V_c^(n+1)
     *
     *     =
     *
     *     L I^n/dt
     *     - (1-theta) R I^n
     *     - (1-theta) AS V_c^n
     *     - (1-theta) U_e^n
     *     - theta U_e^(n+1).
     */
    for (size_t i = 0; i < ne; ++i)
    {
        double value =
        0.0;

        for (size_t j = 0; j < ne; ++j)
        {
            value +=
            system->L->data[
                transient_matrix_index(
                    i,
                    j,
                    ne
                )
            ]
            * edge_current_n[j]
            / system->dt;
        }

        value -=
        one_minus_theta
        * system->R->diagonal[i]
        * edge_current_n[i];

        for (size_t j = 0; j < nv; ++j)
        {
            value -=
            one_minus_theta
            * system->AS.data[
                transient_matrix_index(
                    i,
                    j,
                    nv
                )
            ]
            * node_self_voltage_n[j];
        }

        if (edge_voltage_n != NULL)
        {
            value -=
            one_minus_theta
            * edge_voltage_n[i];
        }

        if (edge_voltage_np1 != NULL)
        {
            value -=
            theta
            * edge_voltage_np1[i];
        }

        system->rhs[
            i0 + i
        ] =
        value;
    }


    /*
     * NODE EQUATION:
     *
     *     -A^T I
     *     + F dV_c/dt
     *
     *     =
     *
     *     J.
     *
     *
     * Theta-схема:
     *
     *     -theta A^T I^(n+1)
     *     + F V_c^(n+1)/dt
     *
     *     =
     *
     *     F V_c^n/dt
     *     + (1-theta) A^T I^n
     *     + (1-theta) J^n
     *     + theta J^(n+1).
     */
    for (size_t j = 0; j < nv; ++j)
    {
        double value =
        system->F[j]
        * node_self_voltage_n[j]
        / system->dt;

        for (size_t e = 0; e < ne; ++e)
        {
            value +=
            one_minus_theta
            * system->A->data[
                transient_matrix_index(
                    e,
                    j,
                    nv
                )
            ]
            * edge_current_n[e];
        }

        if (node_source_n != NULL)
        {
            value +=
            one_minus_theta
            * node_source_n[j];
        }

        if (node_source_np1 != NULL)
        {
            value +=
            theta
            * node_source_np1[j];
        }

        system->rhs[
            v0 + j
        ] =
        value;
    }

    return 0;
}


int transient_build_rhs_stage1(
    TransientSystem *system,
    const double *edge_current_n,
    const double *node_self_voltage_n,
    const double *node_source_n,
    const double *node_source_np1,
    const double *edge_voltage_n,
    const double *edge_voltage_np1
)
{
    if (system == NULL ||
        !system->initialized ||
        system->stage != TRANSIENT_STAGE_1)
    {
        return -1;
    }

    return transient_build_rhs_common(
        system,
        edge_current_n,
        node_self_voltage_n,
        node_source_n,
        node_source_np1,
        edge_voltage_n,
        edge_voltage_np1
    );
}


int transient_build_rhs_stage2(
    TransientSystem *system,
    const double *edge_current_n,
    const double *node_self_voltage_n,
    double shunt_current_n,
    const double *node_source_n,
    const double *node_source_np1,
    const double *edge_voltage_n,
    const double *edge_voltage_np1
)
{
    if (system == NULL ||
        !system->initialized ||
        system->stage != TRANSIENT_STAGE_2 ||
        system->b_R == NULL ||
        !transient_scheme_is_valid(
            system->scheme
        ) ||
        !isfinite(shunt_current_n))
    {
        return -1;
    }

    if (transient_build_rhs_common(
        system,
        edge_current_n,
        node_self_voltage_n,
        node_source_n,
        node_source_np1,
        edge_voltage_n,
        edge_voltage_np1
    ) != 0)
    {
        return -1;
    }

    const double theta =
    transient_scheme_theta(
        system->scheme
    );

    const double one_minus_theta =
    1.0 - theta;

    size_t nv =
    system->n_nodes;

    size_t v0 =
    transient_voltage_offset(
        system
    );

    /*
     * Из:
     *
     *     -A^T I
     *     + F dV_c/dt
     *     - b_R I_shunt
     *
     *     =
     *
     *     I_gamma
     *
     * previous shunt contribution:
     *
     *     +(1-theta) b_R I_shunt^n.
     */
    for (size_t j = 0; j < nv; ++j)
    {
        system->rhs[
            v0 + j
        ] +=
        one_minus_theta
        * system->b_R[j]
        * shunt_current_n;
    }

    size_t s =
    transient_shunt_index(
        system
    );

    system->rhs[s] =
    0.0;

    return 0;
}



int transient_build_rhs_slot_cells(
    TransientSystem *system,
    const double *edge_current_n,
    const double *node_self_voltage_n,
    const double *slot_inductor_current_n,
    double shunt_current_n,
    const double *node_source_n,
    const double *node_source_np1,
    const double *edge_voltage_n,
    const double *edge_voltage_np1
)
{
    if (system == NULL ||
        !system->initialized ||
        system->stage != TRANSIENT_SLOT_CELLS ||
        slot_inductor_current_n == NULL ||
        system->slot_circuit == NULL ||
        system->slot_mass_matrix == NULL ||
        system->slot_voltage_operators == NULL)
    {
        return -1;
    }

    /*
     * Build the ordinary PEEC edge and node RHS first.
     *
     * The node derivative term is corrected below from:
     *
     *     F V_c^n / dt
     *
     * to:
     *
     *     [F + M_slot] V_c^n / dt.
     */
    if (transient_build_rhs_common(
        system,
        edge_current_n,
        node_self_voltage_n,
        node_source_n,
        node_source_np1,
        edge_voltage_n,
        edge_voltage_np1
    ) != 0)
    {
        return -1;
    }

    const double theta =
    transient_scheme_theta(
        system->scheme
    );

    const double one_minus_theta =
    1.0
    - theta;

    size_t nv =
    system->n_nodes;

    size_t v0 =
    transient_voltage_offset(
        system
    );

    size_t s0 =
    transient_slot_current_offset(
        system
    );

    /*
     * Add:
     *
     *     M_slot V_c^n / dt.
     */
    for (size_t row = 0;
         row < nv;
    ++row)
    {
        double value =
        0.0;

        for (size_t col = 0;
             col < nv;
        ++col)
        {
            value +=
            system->slot_mass_matrix[
                row * nv + col
            ]
            * node_self_voltage_n[col]
            / system->dt;
        }

        system->rhs[
            v0 + row
        ] +=
        value;
    }

    /*
     * Previous slot-inductor contribution to node equation:
     *
     *     +(1-theta) B_slot I_L^n.
     *
     * Because:
     *
     *     -B_slot I_L
     *
     * is on the left side.
     */
    for (size_t m = 0;
         m < system->n_slot_cells;
    ++m)
    {
        const SlotCell *cell =
        &system->slot_circuit->cells[m];

        double value =
        one_minus_theta
        * slot_inductor_current_n[m];

        system->rhs[
            v0 + cell->node_a
        ] -=
        value;

        system->rhs[
            v0 + cell->node_b
        ] +=
        value;

        /*
         * Slot-inductor equation:
         *
         *     L_m dI_m/dt + c_m^T V_c = 0.
         *
         * RHS:
         *
         *     L_m I_m^n/dt
         *     - (1-theta)c_m^T V_c^n.
         */
        const double *c =
        &system->slot_voltage_operators[
            m * nv
        ];

        double branch_rhs =
        cell->inductance
        * slot_inductor_current_n[m]
        / system->dt;

        for (size_t j = 0;
             j < nv;
        ++j)
        {
            branch_rhs -=
            one_minus_theta
            * c[j]
            * node_self_voltage_n[j];
        }

        system->rhs[
            s0 + m
        ] =
        branch_rhs;
    }

    if (system->shunt != NULL)
    {
        if (!isfinite(
            shunt_current_n
        ))
        {
            return -1;
        }

        for (size_t j = 0;
             j < nv;
        ++j)
        {
            system->rhs[
                v0 + j
            ] +=
            one_minus_theta
            * system->b_R[j]
            * shunt_current_n;
        }

        system->rhs[
            transient_shunt_index(
                system
            )
        ] =
        0.0;
    }

    return 0;
}


int transient_extract_slot_cells_solution(
    const TransientSystem *system,
    double *edge_current,
    double *node_self_voltage,
    double *slot_inductor_current,
    double *shunt_current
)
{
    if (system == NULL ||
        !system->initialized ||
        system->stage != TRANSIENT_SLOT_CELLS ||
        system->solution == NULL ||
        edge_current == NULL ||
        node_self_voltage == NULL ||
        slot_inductor_current == NULL)
    {
        return -1;
    }

    memcpy(
        edge_current,
        system->solution,
        system->n_edges
        * sizeof(double)
    );

    memcpy(
        node_self_voltage,
        system->solution
        + system->n_edges,
        system->n_nodes
        * sizeof(double)
    );

    memcpy(
        slot_inductor_current,
        system->solution
        + transient_slot_current_offset(
            system
        ),
        system->n_slot_cells
        * sizeof(double)
    );

    if (system->shunt != NULL)
    {
        if (shunt_current == NULL)
        {
            return -1;
        }

        *shunt_current =
        system->solution[
            transient_shunt_index(
                system
            )
        ];
    }

    return 0;
}


int transient_factorize(
    TransientSystem *system,
    int parallel_threads
)
{
    if (system == NULL ||
        !system->initialized ||
        system->system_matrix == NULL ||
        system->n_total == 0 ||
        parallel_threads < 1)
    {
        return -1;
    }

    lu_free(
        &system->lu
    );

    lu_init(
        &system->lu
    );

    LUOptions options;

    lu_options_default(
        &options,
        parallel_threads
    );

    if (lu_factorize(
        system->system_matrix,
        system->n_total,
        &options,
        &system->lu
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: transient LU factorization failed.\n"
        );

        system->factorized =
        0;

        return -1;
    }

    system->factorized =
    1;

    return 0;
}


int transient_solve_step(
    TransientSystem *system
)
{
    if (system == NULL ||
        !system->initialized ||
        !system->factorized ||
        system->rhs == NULL ||
        system->solution == NULL)
    {
        return -1;
    }

    if (lu_solve(
        &system->lu,
        system->rhs,
        system->solution
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: transient linear solve failed.\n"
        );

        return -1;
    }

    return 0;
}


int transient_extract_stage1_solution(
    const TransientSystem *system,
    double *edge_current,
    double *node_self_voltage
)
{
    if (system == NULL ||
        !system->initialized ||
        system->stage != TRANSIENT_STAGE_1 ||
        system->solution == NULL ||
        edge_current == NULL ||
        node_self_voltage == NULL)
    {
        return -1;
    }

    memcpy(
        edge_current,
        system->solution,
        system->n_edges
        * sizeof(double)
    );

    memcpy(
        node_self_voltage,
        system->solution
        + system->n_edges,
        system->n_nodes
        * sizeof(double)
    );

    return 0;
}


int transient_extract_stage2_solution(
    const TransientSystem *system,
    double *edge_current,
    double *node_self_voltage,
    double *shunt_current
)
{
    if (system == NULL ||
        !system->initialized ||
        system->stage != TRANSIENT_STAGE_2 ||
        system->solution == NULL ||
        edge_current == NULL ||
        node_self_voltage == NULL ||
        shunt_current == NULL)
    {
        return -1;
    }

    memcpy(
        edge_current,
        system->solution,
        system->n_edges
        * sizeof(double)
    );

    memcpy(
        node_self_voltage,
        system->solution
        + system->n_edges,
        system->n_nodes
        * sizeof(double)
    );

    *shunt_current =
    system->solution[
        transient_shunt_index(
            system
        )
    ];

    return 0;
}


int transient_recover_node_state(
    const TransientSystem *system,
    const double *node_self_voltage,
    double *node_potential,
    double *node_charge
)
{
    if (system == NULL ||
        !system->initialized ||
        system->P == NULL ||
        node_self_voltage == NULL)
    {
        return -1;
    }

    if (node_potential == NULL &&
        node_charge == NULL)
    {
        return 0;
    }

    if (node_charge != NULL)
    {
        if (potential_self_voltage_to_charge(
            system->P,
            node_self_voltage,
            node_charge
        ) != 0)
        {
            return -1;
        }
    }

    if (node_potential != NULL)
    {
        if (potential_self_voltage_to_full_potential(
            system->P,
            node_self_voltage,
            node_potential
        ) != 0)
        {
            return -1;
        }
    }

    return 0;
}


void transient_print_info(
    const TransientSystem *system
)
{
    if (system == NULL)
    {
        return;
    }

    printf("\n");
    printf("============================================================\n");
    printf("TRANSIENT PEEC / CES SYSTEM\n");
    printf("============================================================\n");

    if (system->stage == TRANSIENT_STAGE_1)
    {
        printf(
            "Stage              : STAGE 1\n"
        );
    }
    else if (system->stage == TRANSIENT_STAGE_2)
    {
        printf(
            "Stage              : STAGE 2 / SHUNT\n"
        );
    }
    else if (system->stage == TRANSIENT_SLOT_CELLS)
    {
        printf(
            "Stage              : ONE-STAGE / SLOT CELLS\n"
        );
    }
    else
    {
        printf(
            "Stage              : UNKNOWN\n"
        );
    }

    if (system->scheme == TRANSIENT_BACKWARD_EULER)
    {
        printf(
            "Time scheme        : BACKWARD EULER\n"
        );

        printf(
            "Temporal order     : 1\n"
        );
    }
    else if (system->scheme == TRANSIENT_TRAPEZOIDAL)
    {
        printf(
            "Time scheme        : TRAPEZOIDAL RULE\n"
        );

        printf(
            "Temporal order     : 2\n"
        );
    }
    else
    {
        printf(
            "Time scheme        : UNKNOWN\n"
        );
    }

    printf(
        "Node unknown       : V_c (self voltage)\n"
    );

    printf(
        "Electrostatic form : F_j=1/P_jj, phi=S*V_c, AS=A*S\n"
    );

    printf(
        "dt                 : %.9e s\n",
        system->dt
    );

    printf(
        "Edge unknowns      : %zu\n",
        system->n_edges
    );

    printf(
        "Node unknowns      : %zu\n",
        system->n_nodes
    );

    printf(
        "Slot unknowns      : %zu\n",
        system->stage == TRANSIENT_SLOT_CELLS
        ? system->n_slot_cells
        : 0
    );

    printf(
        "Shunt unknowns     : %d\n",
        system->shunt != NULL
        ? 1
        : 0
    );

    printf(
        "System size        : %zu x %zu\n",
        system->n_total,
        system->n_total
    );

    printf(
        "LU factorized      : %s\n",
        system->factorized
        ? "YES"
        : "NO"
    );

    if (system->shunt != NULL)
    {
        printf(
            "Shunt resistance   : %.9e Ohm\n",
            system->shunt->resistance
        );

        printf(
            "Shunt node A       : %zu\n",
            system->shunt->node_a
        );

        printf(
            "Shunt node B       : %zu\n",
            system->shunt->node_b
        );
    }

    printf("============================================================\n");
}
