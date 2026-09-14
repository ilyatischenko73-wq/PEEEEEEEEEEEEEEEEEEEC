#ifndef PEEC_SLOT_CIRCUIT_H
#define PEEC_SLOT_CIRCUIT_H

#include <stddef.h>


/*
 * ============================================================
 * DISTRIBUTED EQUIVALENT SLOT CELLS
 * ============================================================
 *
 * The aperture is not replaced by one global LC branch.
 *
 * It is discretized along its length into local equivalent cells.
 * Cell m connects two opposite PEEC nodes:
 *
 *     node_a <---- L_m || C_m ----> node_b
 *
 * and represents a longitudinal support length:
 *
 *     dl_m.
 *
 * The current orientation is:
 *
 *     node_a -> node_b.
 */
typedef struct
{
    size_t node_a;
    size_t node_b;

    double support_length;

    double inductance;
    double capacitance;

} SlotCell;


/*
 * ============================================================
 * QUASI-STATIC SLOT-LINE MODEL
 * ============================================================
 *
 * The narrow rectangular aperture is approximated by a
 * coplanar-strip transmission line.
 *
 * First compute its characteristic impedance Z0s. Then:
 *
 *     v  = 1 / sqrt(mu * epsilon)
 *
 *     L' = Z0s / v
 *
 *     C' = 1 / (Z0s * v)
 *
 * and for every local cell:
 *
 *     L_m = L' * dl_m
 *
 *     C_m = C' * dl_m.
 *
 * slot_width:
 *
 *     physical aperture width w [m].
 *
 * wall_span:
 *
 *     transverse wall dimension b [m] used by the
 *     Robinson/Gupta coplanar-strip approximation.
 *
 * wall_thickness:
 *
 *     conducting wall thickness t [m].
 */
typedef struct
{
    double slot_width;
    double wall_span;
    double wall_thickness;

    double epsilon;
    double mu;

    double effective_width;

    double characteristic_impedance;
    double wave_velocity;

    double inductance_per_length;
    double capacitance_per_length;

} SlotLineModel;


typedef struct
{
    size_t n_cells;
    SlotCell *cells;

    SlotLineModel line;

} SlotCircuit;


void slot_circuit_init(
    SlotCircuit *circuit
);

void slot_circuit_free(
    SlotCircuit *circuit
);


/*
 * Build the geometry-based quasi-static approximation.
 */
int slot_line_model_build(
    SlotLineModel *model,
    double slot_width,
    double wall_span,
    double wall_thickness,
    double epsilon,
    double mu
);


/*
 * Text file format:
 *
 *     # node_a node_b support_length_m
 *     179 187 0.25
 *     180 186 0.25
 *
 * Node numbers are INTERNAL 0-based PEEC node indices.
 */
int slot_circuit_load_cells(
    SlotCircuit *circuit,
    const char *filename,
    const SlotLineModel *model,
    size_t n_mesh_nodes
);


void slot_circuit_print_info(
    const SlotCircuit *circuit
);

#endif
