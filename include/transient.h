#ifndef PEEC_TRANSIENT_H
#define PEEC_TRANSIENT_H

#include "incidence.h"
#include "inductance.h"
#include "internal_circuit.h"
#include "lu.h"
#include "matrix.h"
#include "potential.h"
#include "resistance.h"
#include "slot_circuit.h"

#include <stddef.h>


/*
 * ============================================================
 * РЕЖИМ ВРЕМЕННОЙ ЗАДАЧИ
 * ============================================================
 */
typedef enum
{
    TRANSIENT_STAGE_1 = 0,
    TRANSIENT_STAGE_2 = 1,
    TRANSIENT_SLOT_CELLS = 2

} TransientStage;


/*
 * ============================================================
 * ВРЕМЕННАЯ СХЕМА
 * ============================================================
 *
 * BACKWARD EULER:
 *
 *     theta = 1,
 *
 *     первый порядок.
 *
 *
 * TRAPEZOIDAL:
 *
 *     theta = 1/2,
 *
 *     второй порядок.
 */
typedef enum
{
    TRANSIENT_BACKWARD_EULER = 1,
    TRANSIENT_TRAPEZOIDAL = 2

} TransientScheme;


/*
 * ============================================================
 * TRANSIENT PEEC / ЧЭС SYSTEM
 * ============================================================
 *
 * Неизвестные:
 *
 *     I   - токи глобальных ребер;
 *
 *     V_c - собственные узловые напряжения.
 *
 *
 * Электростатические операторы:
 *
 *     F_j = 1 / P_jj,
 *
 *     S_ja = P_ja / P_aa,
 *
 *     phi = S V_c,
 *
 *     q   = F V_c.
 *
 *
 * STAGE 1:
 *
 *     L dI/dt + R I + A S V_c + U_e = 0,
 *
 *     A^T I - F dV_c/dt + J = 0.
 *
 *
 * После умножения второго уравнения на -1:
 *
 *     -A^T I + F dV_c/dt = J.
 *
 *
 * STAGE 2:
 *
 *     L dI/dt + R I + A S V_c + U_e = 0,
 *
 *     -A^T I
 *     + F dV_c/dt
 *     - b_R I_shunt
 *
 *     =
 *
 *     I_gamma,
 *
 *
 *     b_R^T phi + R_shunt I_shunt = 0,
 *
 * где:
 *
 *     phi = S V_c.
 *
 *
 * Таким образом алгебраическое уравнение шунта в переменной
 * V_c использует:
 *
 *     b_R^T S V_c + R_shunt I_shunt = 0.
 */
typedef struct
{
    TransientStage stage;
    TransientScheme scheme;

    size_t n_edges;
    size_t n_nodes;
    size_t n_total;

    double dt;

    /*
     * Внешние PEEC-матрицы.
     *
     * Память принадлежит вызывающей стороне.
     */
    const InductanceMatrix *L;
    const ResistanceMatrix *R;
    const IncidenceMatrix *A;
    const PotentialMatrix *P;

    /*
     * --------------------------------------------------------
     * ЭЛЕКТРОСТАТИЧЕСКИЕ ОПЕРАТОРЫ
     * --------------------------------------------------------
     *
     * F:
     *
     *     диагональ F_j = 1/P_jj.
     *
     * Храним только диагональ.
     */
    double *F;

    /*
     * AS:
     *
     *     Ne x Nv.
     *
     *     AS = A S,
     *
     *     S_ja = P_ja/P_aa.
     */
    Matrix AS;

    /*
     * --------------------------------------------------------
     * STAGE 2
     * --------------------------------------------------------
     */
    const InternalCircuitShunt *shunt;

    /*
     * b_R:
     *
     *     b_R[A] = -1,
     *     b_R[B] = +1.
     */
    double *b_R;

    /*
     * shunt_voltage_operator:
     *
     *     c^T = b_R^T S.
     *
     *
     * Тогда:
     *
     *     V_B(full) - V_A(full)
     *
     *     =
     *
     *     c^T V_c.
     */
    double *shunt_voltage_operator;

    /*
     * --------------------------------------------------------
     * ONE-STAGE DISTRIBUTED SLOT CELLS
     * --------------------------------------------------------
     */
    const SlotCircuit *slot_circuit;

    size_t n_slot_cells;

    /*
     * Dense Nv x Nv contribution:
     *
     *     M_slot = sum_m C_m b_m c_m^T
     *
     * so the node derivative operator is:
     *
     *     diag(F) + M_slot.
     */
    double *slot_mass_matrix;

    /*
     * Row m:
     *
     *     c_m^T = b_m^T S
     *
     * size:
     *
     *     Nslot x Nv.
     */
    double *slot_voltage_operators;

    /*
     * --------------------------------------------------------
     * ЛИНЕЙНАЯ СИСТЕМА
     * --------------------------------------------------------
     */
    double *system_matrix;
    double *rhs;
    double *solution;

    LUFactorization lu;

    int factorized;
    int initialized;

} TransientSystem;


void transient_system_init(
    TransientSystem *system
);

void transient_system_free(
    TransientSystem *system
);


/*
 * ============================================================
 * СОЗДАНИЕ STAGE 1
 * ============================================================
 *
 * Общая theta-матрица:
 *
 *     [ L/dt + theta R       theta A S ]
 *     [                                  ]
 *     [ -theta A^T           F/dt       ].
 */
int transient_system_create_stage1(
    TransientSystem *system,
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const PotentialMatrix *P,
    TransientScheme scheme,
    double dt
);


/*
 * ============================================================
 * СОЗДАНИЕ STAGE 2
 * ============================================================
 *
 *     [ L/dt + theta R   theta AS     0            ]
 *     [                                             ]
 *     [ -theta A^T       F/dt        -theta b_R    ]
 *     [                                             ]
 *     [ 0                b_R^T S      R_shunt      ].
 */
int transient_system_create_stage2(
    TransientSystem *system,
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const PotentialMatrix *P,
    const InternalCircuitShunt *shunt,
    TransientScheme scheme,
    double dt
);


/*
 * Строит:
 *
 *     b_R
 *
 * и:
 *
 *     b_R^T S.
 */
int transient_build_shunt_operators(
    TransientSystem *system
);


int transient_build_system_matrix(
    TransientSystem *system
);


/*
 * ============================================================
 * ONE-STAGE SLOT-CELL SYSTEM
 * ============================================================
 *
 * Unknowns:
 *
 *     [ I_edges, V_c, I_L,slot, (optional I_shunt) ].
 *
 * The shunt pointer may be NULL.
 */
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
);


size_t transient_slot_current_offset(
    const TransientSystem *system
);


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
);


int transient_extract_slot_cells_solution(
    const TransientSystem *system,
    double *edge_current,
    double *node_self_voltage,
    double *slot_inductor_current,
    double *shunt_current
);


int transient_factorize(
    TransientSystem *system,
    int parallel_threads
);


/*
 * ============================================================
 * RHS STAGE 1
 * ============================================================
 *
 * node_source:
 *
 *     J(t) [A].
 *
 *
 * Для lightning-current:
 *
 *     J = I_s.
 *
 *
 * Для incident-pulse:
 *
 *     J = 0,
 *
 *     U_e != 0.
 */
int transient_build_rhs_stage1(
    TransientSystem *system,
    const double *edge_current_n,
    const double *node_self_voltage_n,
    const double *node_source_n,
    const double *node_source_np1,
    const double *edge_voltage_n,
    const double *edge_voltage_np1
);


/*
 * ============================================================
 * RHS STAGE 2
 * ============================================================
 *
 * node_source:
 *
 *     I_gamma.
 */
int transient_build_rhs_stage2(
    TransientSystem *system,
    const double *edge_current_n,
    const double *node_self_voltage_n,
    double shunt_current_n,
    const double *node_source_n,
    const double *node_source_np1,
    const double *edge_voltage_n,
    const double *edge_voltage_np1
);


int transient_solve_step(
    TransientSystem *system
);


size_t transient_current_offset(
    const TransientSystem *system
);

size_t transient_voltage_offset(
    const TransientSystem *system
);

size_t transient_shunt_index(
    const TransientSystem *system
);


int transient_extract_stage1_solution(
    const TransientSystem *system,
    double *edge_current,
    double *node_self_voltage
);


int transient_extract_stage2_solution(
    const TransientSystem *system,
    double *edge_current,
    double *node_self_voltage,
    double *shunt_current
);


/*
 * ============================================================
 * ВОССТАНОВЛЕНИЕ ФИЗИЧЕСКИХ УЗЛОВЫХ ВЕЛИЧИН
 * ============================================================
 *
 * Из внутренней неизвестной:
 *
 *     V_c
 *
 * получаем:
 *
 *     q   = F V_c,
 *
 *     phi = S V_c = P q.
 *
 *
 * node_potential и node_charge могут быть NULL независимо.
 */
int transient_recover_node_state(
    const TransientSystem *system,
    const double *node_self_voltage,
    double *node_potential,
    double *node_charge
);


void transient_print_info(
    const TransientSystem *system
);


#endif
