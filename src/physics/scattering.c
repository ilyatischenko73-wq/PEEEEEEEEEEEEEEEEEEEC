#include "scattering.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>


/*
 * ============================================================
 * sinc(x) = sin(x) / x
 * ============================================================
 *
 * При x -> 0 прямое вычисление sin(x)/x приводит
 * к потере точности.
 *
 * Поэтому около нуля используем ряд:
 *
 *     sinc(x)
 *
 *     = 1 - x^2/6 + x^4/120 - ...
 */
static double scattering_sinc(double x)
{
    double ax = fabs(x);

    if (ax < 1.0e-4)
    {
        double x2 = x * x;
        double x4 = x2 * x2;

        return 1.0 - x2 / 6.0 + x4 / 120.0;
    }

    return sin(x) / x;
}


/*
 * ============================================================
 * ПАРАМЕТРЫ ПО УМОЛЧАНИЮ
 * ============================================================
 */
void scattering_options_default(
    ScatteringOptions *options,
    int parallel_threads
)
{
    if (options == NULL)
    {
        return;
    }

    options->integration =
        SCATTERING_EDGE_MIDPOINT;

    options->parallel_threads =
        parallel_threads;
}


/*
 * ============================================================
 * ПРОВЕРКА СЕТКИ
 * ============================================================
 */
static int scattering_check_mesh(
    const Mesh *mesh
)
{
    if (mesh == NULL)
    {
        return -1;
    }

    if (mesh->n_edges == 0)
    {
        fprintf(
            stderr,
            "ERROR: scattering mesh contains no edges.\n"
        );

        return -1;
    }

    if (mesh->edge_vectors == NULL ||
        mesh->edge_centers == NULL)
    {
        fprintf(
            stderr,
            "ERROR: scattering requires edge vectors and edge centers.\n"
        );

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ПРОВЕРКА ПАРАМЕТРОВ
 * ============================================================
 */
static int scattering_check_options(
    const ScatteringOptions *options
)
{
    if (options == NULL)
    {
        return -1;
    }

    if (options->parallel_threads < 1)
    {
        fprintf(
            stderr,
            "ERROR: invalid OpenMP thread count in scattering.\n"
        );

        return -1;
    }

    if (options->integration != SCATTERING_EDGE_MIDPOINT &&
        options->integration != SCATTERING_EDGE_ANALYTIC_PLANE_WAVE)
    {
        fprintf(
            stderr,
            "ERROR: unsupported scattering edge integration method.\n"
        );

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ВЕЩЕСТВЕННОЕ СКАЛЯРНОЕ ПРОИЗВЕДЕНИЕ
 * ============================================================
 */
static inline double dot_real(
    const double a[3],
    const double b[3]
)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


/*
 * ============================================================
 * КОМПЛЕКСНЫЙ ВЕКТОР . ВЕЩЕСТВЕННЫЙ ВЕКТОР
 * ============================================================
 */
static inline Complex dot_complex_real(
    const Complex a[3],
    const double b[3]
)
{
    Complex value;

    value.re = a[0].re * b[0] + a[1].re * b[1] + a[2].re * b[2];
    value.im = a[0].im * b[0] + a[1].im * b[1] + a[2].im * b[2];

    return value;
}


/*
 * ============================================================
 * MIDPOINT ДЛЯ ГАРМОНИЧЕСКОГО ПОЛЯ
 * ============================================================
 *
 *     U_e ~= E_hat(r_c) . l_e.
 */
static int scattering_harmonic_midpoint(
    const IncidentField *field,
    const double center[3],
    const double edge_vector[3],
    Complex *value
)
{
    Complex E[3];

    if (incident_field_evaluate_phasor(
            field,
            center,
            E
        ) != 0)
    {
        return -1;
    }

    *value = dot_complex_real(
        E,
        edge_vector
    );

    return 0;
}


/*
 * ============================================================
 * ТОЧНЫЙ ИНТЕГРАЛ ПЛОСКОЙ ВОЛНЫ ПО РЕБРУ
 * ============================================================
 *
 * Плоская гармоническая волна:
 *
 *     E_hat(r) =
 *
 *     E0 exp(
 *         -j*k0*k_hat.r
 *         +j*phase
 *     ).
 *
 *
 * Ребро:
 *
 *     r(s) = r_c + s*l,
 *
 *     -1/2 <= s <= 1/2.
 *
 *
 * Тогда:
 *
 *     dl = l ds
 *
 * и:
 *
 *     U =
 *
 *     integral E_hat(r(s)).l ds.
 *
 *
 * После интегрирования:
 *
 *     U =
 *
 *     [E_hat(r_c).l]
 *
 *     * sinc(q/2),
 *
 * где:
 *
 *     q = k0 * k_hat.l.
 *
 *
 * Если:
 *
 *     k_hat.l = 0,
 *
 * то:
 *
 *     sinc(0) = 1
 *
 * и точный результат совпадает с midpoint.
 */
static int scattering_harmonic_analytic_plane_wave(
    const IncidentField *field,
    const double center[3],
    const double edge_vector[3],
    Complex *value
)
{
    Complex E_center[3];

    if (incident_field_evaluate_phasor(
            field,
            center,
            E_center
        ) != 0)
    {
        return -1;
    }

    /*
     * Поле в центре ребра, спроецированное
     * на ориентированный вектор ребра.
     *
     * Размерность:
     *
     *     [V/m] * [m] = [V].
     */
    Complex midpoint =
        dot_complex_real(
            E_center,
            edge_vector
        );


    /*
     * Фазовое изменение вдоль всего ребра:
     *
     *     q = k0 * k_hat.l.
     *
     * q безразмерен [rad].
     */
    double k_dot_l =
        dot_real(
            field->k_hat,
            edge_vector
        );

    double q =
        field->wave_number
        * k_dot_l;


    /*
     * Точный поправочный множитель.
     */
    double factor =
        scattering_sinc(
            0.5 * q
        );


    value->re =
        midpoint.re * factor;

    value->im =
        midpoint.im * factor;

    return 0;
}


/*
 * ============================================================
 * ГАРМОНИЧЕСКОЕ ВОЗБУЖДЕНИЕ
 * ============================================================
 */
int scattering_build_harmonic_excitation(
    const Mesh *mesh,
    const IncidentField *field,
    const ScatteringOptions *options,
    Complex *edge_voltage
)
{
    if (mesh == NULL ||
        field == NULL ||
        options == NULL ||
        edge_voltage == NULL)
    {
        return -1;
    }

    if (scattering_check_mesh(mesh) != 0 ||
        scattering_check_options(options) != 0)
    {
        return -1;
    }


    /*
     * Импульс не является одночастотным фазором.
     */
    if (field->wave_type == INCIDENT_WAVE_PULSE)
    {
        fprintf(
            stderr,
            "ERROR: pulse incident field cannot be used directly "
            "in harmonic scattering.\n"
        );

        return -1;
    }


    size_t ne =
        mesh->n_edges;


    /*
     * Внутри parallel-for не делаем fprintf.
     *
     * При ошибке записываем NaN, после чего
     * проверяем результат последовательно.
     */
    #pragma omp parallel for schedule(static) num_threads(options->parallel_threads)
    for (long long e_signed = 0;
         e_signed < (long long)ne;
         ++e_signed)
    {
        size_t e =
            (size_t)e_signed;

        const double *center =
            &mesh->edge_centers[3 * e];

        const double *edge_vector =
            &mesh->edge_vectors[3 * e];

        Complex value;

        int status = -1;


        /*
         * --------------------------------------------
         * MIDPOINT
         * --------------------------------------------
         */
        if (options->integration ==
            SCATTERING_EDGE_MIDPOINT)
        {
            status =
                scattering_harmonic_midpoint(
                    field,
                    center,
                    edge_vector,
                    &value
                );
        }


        /*
         * --------------------------------------------
         * ТОЧНЫЙ ИНТЕГРАЛ ПЛОСКОЙ ВОЛНЫ
         * --------------------------------------------
         */
        else if (options->integration ==
                 SCATTERING_EDGE_ANALYTIC_PLANE_WAVE)
        {
            status =
                scattering_harmonic_analytic_plane_wave(
                    field,
                    center,
                    edge_vector,
                    &value
                );
        }


        if (status != 0)
        {
            edge_voltage[e].re = NAN;
            edge_voltage[e].im = NAN;

            continue;
        }

        edge_voltage[e] =
            value;
    }


    for (size_t e = 0; e < ne; ++e)
    {
        if (!isfinite(edge_voltage[e].re) ||
            !isfinite(edge_voltage[e].im))
        {
            fprintf(
                stderr,
                "ERROR: invalid harmonic incident excitation "
                "at edge %zu.\n",
                e
            );

            return -1;
        }
    }

    return 0;
}


/*
 * ============================================================
 * ВРЕМЕННОЕ ВОЗБУЖДЕНИЕ
 * ============================================================
 *
 * Для временной задачи пока всегда используем midpoint:
 *
 *     U_e(t) ~= E(r_c,t).l_e.
 *
 *
 * Причина:
 *
 * аналитический sinc выше получен именно для
 * монохроматического комплексного фазора.
 *
 * Для COS/SIN временная реализация тоже могла бы иметь
 * аналитический интеграл, но для PULSE уже потребуется
 * интегрировать задержанное поле вдоль ребра.
 */
int scattering_build_time_excitation(
    const Mesh *mesh,
    const IncidentField *field,
    const ScatteringOptions *options,
    double time,
    double *edge_voltage
)
{
    if (mesh == NULL ||
        field == NULL ||
        options == NULL ||
        edge_voltage == NULL ||
        !isfinite(time))
    {
        return -1;
    }

    if (scattering_check_mesh(mesh) != 0 ||
        scattering_check_options(options) != 0)
    {
        return -1;
    }

    size_t ne =
        mesh->n_edges;


    #pragma omp parallel for schedule(static) num_threads(options->parallel_threads)
    for (long long e_signed = 0;
         e_signed < (long long)ne;
         ++e_signed)
    {
        size_t e =
            (size_t)e_signed;

        const double *center =
            &mesh->edge_centers[3 * e];

        const double *edge_vector =
            &mesh->edge_vectors[3 * e];

        double E[3];

        int status =
            incident_field_evaluate_time(
                field,
                center,
                time,
                E
            );

        if (status != 0)
        {
            edge_voltage[e] = NAN;
            continue;
        }


        /*
         *     U_e(t)
         *
         *     ~= E(r_c,t).l_e.
         */
        edge_voltage[e] =
            dot_real(
                E,
                edge_vector
            );
    }


    for (size_t e = 0; e < ne; ++e)
    {
        if (!isfinite(edge_voltage[e]))
        {
            fprintf(
                stderr,
                "ERROR: invalid time-domain incident excitation "
                "at edge %zu.\n",
                e
            );

            return -1;
        }
    }

    return 0;
}


/*
 * ============================================================
 * НУЛЕВОЙ УЗЛОВОЙ ТОК
 * ============================================================
 */
void scattering_zero_node_current_real(
    size_t n_nodes,
    double *node_current
)
{
    if (node_current == NULL)
    {
        return;
    }

    memset(
        node_current,
        0,
        n_nodes * sizeof(double)
    );
}


void scattering_zero_node_current_complex(
    size_t n_nodes,
    Complex *node_current
)
{
    if (node_current == NULL)
    {
        return;
    }

    memset(
        node_current,
        0,
        n_nodes * sizeof(Complex)
    );
}


/*
 * ============================================================
 * ДИАГНОСТИКА ГАРМОНИЧЕСКОГО ВОЗБУЖДЕНИЯ
 * ============================================================
 */
void scattering_print_harmonic_excitation_info(
    const Complex *edge_voltage,
    size_t n_edges
)
{
    if (edge_voltage == NULL ||
        n_edges == 0)
    {
        return;
    }

    double min_abs =
        hypot(
            edge_voltage[0].re,
            edge_voltage[0].im
        );

    double max_abs =
        min_abs;

    double sum_squared =
        0.0;


    for (size_t e = 0; e < n_edges; ++e)
    {
        double magnitude =
            hypot(
                edge_voltage[e].re,
                edge_voltage[e].im
            );

        if (magnitude < min_abs)
        {
            min_abs = magnitude;
        }

        if (magnitude > max_abs)
        {
            max_abs = magnitude;
        }

        sum_squared +=
            magnitude * magnitude;
    }


    double rms =
        sqrt(
            sum_squared / (double)n_edges
        );


    printf("\n");
    printf("========================================\n");
    printf("HARMONIC INCIDENT EXCITATION\n");
    printf("========================================\n");

    printf(
        "Edges             : %zu\n",
        n_edges
    );

    printf(
        "Minimum |U_inc|   : %.9e V\n",
        min_abs
    );

    printf(
        "Maximum |U_inc|   : %.9e V\n",
        max_abs
    );

    printf(
        "RMS |U_inc|       : %.9e V\n",
        rms
    );

    printf("========================================\n");
}


/*
 * ============================================================
 * ДИАГНОСТИКА ВРЕМЕННОГО ВОЗБУЖДЕНИЯ
 * ============================================================
 */
void scattering_print_time_excitation_info(
    const double *edge_voltage,
    size_t n_edges,
    double time
)
{
    if (edge_voltage == NULL ||
        n_edges == 0)
    {
        return;
    }

    double min_value =
        edge_voltage[0];

    double max_value =
        edge_voltage[0];

    double max_abs =
        fabs(
            edge_voltage[0]
        );

    double sum_squared =
        0.0;


    for (size_t e = 0; e < n_edges; ++e)
    {
        double value =
            edge_voltage[e];

        double magnitude =
            fabs(value);

        if (value < min_value)
        {
            min_value = value;
        }

        if (value > max_value)
        {
            max_value = value;
        }

        if (magnitude > max_abs)
        {
            max_abs = magnitude;
        }

        sum_squared +=
            value * value;
    }


    double rms =
        sqrt(
            sum_squared / (double)n_edges
        );


    printf("\n");
    printf("========================================\n");
    printf("TIME-DOMAIN INCIDENT EXCITATION\n");
    printf("========================================\n");

    printf(
        "Time              : %.9e s\n",
        time
    );

    printf(
        "Edges             : %zu\n",
        n_edges
    );

    printf(
        "Minimum U_inc     : %.9e V\n",
        min_value
    );

    printf(
        "Maximum U_inc     : %.9e V\n",
        max_value
    );

    printf(
        "Maximum |U_inc|   : %.9e V\n",
        max_abs
    );

    printf(
        "RMS U_inc         : %.9e V\n",
        rms
    );

    printf("========================================\n");
}