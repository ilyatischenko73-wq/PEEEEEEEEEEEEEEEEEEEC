#ifndef PEEC_HARMONIC_H
#define PEEC_HARMONIC_H

#include "complex_lu.h"
#include "complex_matrix.h"
#include "incidence.h"
#include "inductance.h"
#include "matrix.h"
#include "potential.h"
#include "resistance.h"

#include <stddef.h>


/*
 * ============================================================
 * ГАРМОНИЧЕСКАЯ PEEC-СИСТЕМА
 * ============================================================
 *
 * Используем соглашение:
 *
 *     exp(j*omega*t).
 *
 *
 * Классическая PEEC-система:
 *
 *     (R + j*omega*L) I + A V = -U_e
 *
 *     -P A^T I + j*omega V = P I_s.
 *
 *
 * В блочной форме:
 *
 *     [ R+j*w*L      A     ] [ I ]   [ -U_e  ]
 *     [                    ] [   ] = [       ]
 *     [ -P*A^T      j*w*Iv ] [ V ]   [ P I_s ].
 *
 *
 * Неизвестные:
 *
 *     I - токи ребер [A], размер Ne;
 *
 *     V - узловые потенциалы [V], размер Nv.
 *
 *
 * Полный размер:
 *
 *     N = Ne + Nv.
 */
typedef struct
{
    size_t n_edges;
    size_t n_nodes;
    size_t system_size;

    /*
     * Частота [Hz].
     */
    double frequency;

    /*
     * Круговая частота:
     *
     *     omega = 2*pi*f [rad/s].
     */
    double omega;

    /*
     * Полная комплексная матрица:
     *
     *     Z : N x N.
     */
    ComplexMatrix matrix;

    /*
     * LU-факторизация полной системы.
     */
    ComplexLUFactorization lu;

    int matrix_built;
    int factorized;

} HarmonicSystem;


/*
 * ============================================================
 * ГАРМОНИЧЕСКОЕ ВОЗБУЖДЕНИЕ
 * ============================================================
 *
 * edge_voltage:
 *
 *     U_e [V], размер Ne.
 *
 *
 * node_current:
 *
 *     I_s [A], размер Nv.
 *
 *
 * Для обычного рассеяния:
 *
 *     node_current = NULL
 *
 * означает:
 *
 *     I_s = 0.
 *
 *
 * Если:
 *
 *     edge_voltage = NULL,
 *
 * то считается:
 *
 *     U_e = 0.
 */
typedef struct
{
    size_t n_edges;
    size_t n_nodes;

    const Complex *edge_voltage;
    const Complex *node_current;

} HarmonicExcitation;


/*
 * ============================================================
 * ГАРМОНИЧЕСКОЕ РЕШЕНИЕ
 * ============================================================
 */
typedef struct
{
    size_t n_edges;
    size_t n_nodes;

    /*
     * I_e [A].
     */
    Complex *edge_current;

    /*
     * V_j [V].
     */
    Complex *node_voltage;

} HarmonicSolution;


/*
 * ============================================================
 * ДИАГНОСТИКА НЕВЯЗКИ
 * ============================================================
 *
 * Проверяем:
 *
 *     Z x = b,
 *
 * где:
 *
 *     x = [I ; V].
 *
 *
 * ВАЖНО:
 *
 * Верхний и нижний блоки PEEC-системы имеют разные
 * физические размерности:
 *
 * верхний блок:
 *
 *     (R+j*w*L)I + AV = -U
 *
 * имеет размерность [V],
 *
 * а нижний:
 *
 *     -P A^T I + j*w V = P I_s
 *
 * имеет размерность [V/s].
 *
 *
 * Поэтому одна абсолютная невязка всей системы является
 * только вспомогательной диагностикой.
 *
 *
 * Основная величина:
 *
 *     backward_error =
 *
 *                    |r_i|
 *     max_i -----------------------------
 *           sum_j |Z_ij| |x_j| + |b_i|.
 *
 *
 * Это компонентно нормированная backward error.
 */
typedef struct
{
    /*
     * Максимальная абсолютная невязка среди всех строк.
     */
    double absolute_residual;

    /*
     * Старое отношение:
     *
     *     ||Zx-b||_inf / ||b||_inf.
     *
     * Оставляется как дополнительная диагностика.
     */
    double rhs_relative_residual;

    /*
     * Основная нормированная backward error.
     */
    double backward_error;

    /*
     * Отдельная абсолютная невязка верхнего PEEC-блока:
     *
     *     (R+j*w*L)I + AV = -U_e.
     */
    double edge_absolute_residual;

    /*
     * Отдельная абсолютная невязка нижнего PEEC-блока:
     *
     *     -P A^T I + j*w V = P I_s.
     */
    double node_absolute_residual;

    /*
     * Нормированная backward error верхнего блока.
     */
    double edge_backward_error;

    /*
     * Нормированная backward error нижнего блока.
     */
    double node_backward_error;

} HarmonicResidual;


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ / ОСВОБОЖДЕНИЕ
 * ============================================================
 */
void harmonic_system_init(
    HarmonicSystem *system
);


void harmonic_system_free(
    HarmonicSystem *system
);


void harmonic_solution_init(
    HarmonicSolution *solution
);


void harmonic_solution_free(
    HarmonicSolution *solution
);


/*
 * ============================================================
 * ПОСТРОЕНИЕ ПОЛНОЙ PEEC-МАТРИЦЫ
 * ============================================================
 *
 * Строит:
 *
 *     [ R+j*w*L      A     ]
 *     [                    ]
 *     [ -P*A^T      j*w*Iv ].
 *
 *
 * PAT должен быть заранее построен:
 *
 *     PAT = P A^T.
 */
int harmonic_build_system(
    const InductanceMatrix *L,
    const ResistanceMatrix *R,
    const IncidenceMatrix *A,
    const Matrix *PAT,
    double frequency,
    int parallel_threads,
    HarmonicSystem *system
);


/*
 * ============================================================
 * LU-ФАКТОРИЗАЦИЯ
 * ============================================================
 */
int harmonic_factorize(
    HarmonicSystem *system,
    int parallel_threads
);


/*
 * ============================================================
 * ПРАВАЯ ЧАСТЬ
 * ============================================================
 *
 * Строит:
 *
 *     rhs =
 *
 *     [ -U_e  ]
 *     [ P I_s ].
 *
 *
 * Если edge_voltage == NULL:
 *
 *     U_e = 0.
 *
 * Если node_current == NULL:
 *
 *     I_s = 0.
 */
int harmonic_build_rhs(
    const PotentialMatrix *P,
    const HarmonicExcitation *excitation,
    Complex *rhs,
    int parallel_threads
);


/*
 * ============================================================
 * РЕШЕНИЕ
 * ============================================================
 */
int harmonic_solve(
    const HarmonicSystem *system,
    const PotentialMatrix *P,
    const HarmonicExcitation *excitation,
    HarmonicSolution *solution,
    int parallel_threads
);


/*
 * ============================================================
 * ПРОВЕРКА НЕВЯЗКИ
 * ============================================================
 */
int harmonic_compute_residual(
    const HarmonicSystem *system,
    const PotentialMatrix *P,
    const HarmonicExcitation *excitation,
    const HarmonicSolution *solution,
    int parallel_threads,
    HarmonicResidual *residual
);


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void harmonic_print_info(
    const HarmonicSystem *system
);


#endif