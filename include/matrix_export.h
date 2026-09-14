#ifndef PEEC_MATRIX_EXPORT_H
#define PEEC_MATRIX_EXPORT_H

#include "mesh.h"
#include "inductance.h"
#include "potential.h"
#include "resistance.h"
#include "matrix.h"

#include <stddef.h>

/*
 * Сохраняем PEEC-матрицы в Matrix Market.
 *
 * Формат хорошо читается:
 *
 *     Python / scipy.io.mmread
 *     MATLAB
 *     GNU Octave
 *
 * L, P, PAT сохраняются как dense "array".
 * A и R сохраняются как sparse "coordinate".
 */

int matrix_export_dense_mtx(
    const char *filename,
    const double *data,
    size_t rows,
    size_t cols
);

int matrix_export_incidence_mtx(
    const char *filename,
    const Mesh *mesh
);

int matrix_export_resistance_mtx(
    const char *filename,
    const ResistanceMatrix *R
);

/*
 * Сохраняет комплект:
 *
 *     L.mtx
 *     A.mtx
 *     P.mtx
 *     PAT.mtx
 *     R.mtx
 *
 * в directory.
 */
int matrix_export_peec_set(
    const char *directory,
    const Mesh *mesh,
    const InductanceMatrix *L,
    const PotentialMatrix *P,
    const Matrix *PAT,
    const ResistanceMatrix *R
);

#endif
