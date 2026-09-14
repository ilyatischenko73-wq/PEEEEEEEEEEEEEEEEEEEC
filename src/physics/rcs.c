#include "rcs.h"

#include "quadrature.h"

#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


/*
 * ============================================================
 * ФИЗИЧЕСКИЕ КОНСТАНТЫ
 * ============================================================
 */

static const double RCS_MU0 =
    1.25663706212e-6;

static const double RCS_EPS0 =
    8.8541878128e-12;


/*
 * ============================================================
 * ВЕКТОРНАЯ ГЕОМЕТРИЯ
 * ============================================================
 */

static inline double rcs_dot3(
    const double a[3],
    const double b[3]
)
{
    return a[0] * b[0]
         + a[1] * b[1]
         + a[2] * b[2];
}


static inline double rcs_norm3(
    const double a[3]
)
{
    return sqrt(
        a[0] * a[0]
        + a[1] * a[1]
        + a[2] * a[2]
    );
}


static void rcs_cross3(
    const double a[3],
    const double b[3],
    double c[3]
)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}


/*
 * ============================================================
 * КОМПЛЕКСНЫЕ ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
 * ============================================================
 */

static inline Complex rcs_complex_mul(
    Complex a,
    Complex b
)
{
    Complex value;

    value.re = a.re * b.re - a.im * b.im;
    value.im = a.re * b.im + a.im * b.re;

    return value;
}


/*
 *     exp(j*phase).
 */
static inline Complex rcs_complex_exp_j(
    double phase
)
{
    Complex value;

    value.re = cos(phase);
    value.im = sin(phase);

    return value;
}


/*
 * ============================================================
 * QUADPATCH -> POINT
 * ============================================================
 *
 * Та же структура хранения, что используется
 * в dual_mesh.h:
 *
 *     x0 y0 z0 x1 y1 z1 ...
 */
static void rcs_patch_get_point(
    const QuadPatch *patch,
    int index,
    double point[3]
)
{
    point[0] = patch->points[3 * index + 0];
    point[1] = patch->points[3 * index + 1];
    point[2] = patch->points[3 * index + 2];
}


/*
 * ============================================================
 * ПАРАМЕТРИЗАЦИЯ ТРЕУГОЛЬНИКА
 * ============================================================
 *
 * Это та же параметризация, что используется
 * в quadrature.c:
 *
 *     u,v in [0,1],
 *
 *     xi  = u,
 *     eta = (1-u)v.
 *
 *
 *     r =
 *
 *     p0
 *     + xi  (p1-p0)
 *     + eta (p2-p0).
 *
 *
 * Поверхностный Jacobian:
 *
 *     J =
 *
 *     |(p1-p0)x(p2-p0)|
 *     * (1-u).
 */
static void rcs_triangle_map(
    const double p0[3],
    const double p1[3],
    const double p2[3],
    double u,
    double v,
    double r[3],
    double *jacobian
)
{
    double a[3] = {
        p1[0] - p0[0],
        p1[1] - p0[1],
        p1[2] - p0[2]
    };

    double b[3] = {
        p2[0] - p0[0],
        p2[1] - p0[1],
        p2[2] - p0[2]
    };

    double xi =
        u;

    double eta =
        (1.0 - u) * v;


    r[0] =
        p0[0]
        + xi * a[0]
        + eta * b[0];

    r[1] =
        p0[1]
        + xi * a[1]
        + eta * b[1];

    r[2] =
        p0[2]
        + xi * a[2]
        + eta * b[2];


    double cross[3];

    rcs_cross3(
        a,
        b,
        cross
    );

    *jacobian =
        rcs_norm3(cross)
        * (1.0 - u);
}


/*
 * ============================================================
 * ИНТЕГРАЛ ПО ТРЕУГОЛЬНИКУ
 * ============================================================
 *
 * Вычисляет:
 *
 *     G_T =
 *
 *     integral_T
 *
 *     exp(
 *         j*k*r_hat.r
 *     )
 *
 *     dS.
 *
 *
 * В отличие от интеграла 1/R здесь НЕТ сингулярности.
 *
 * Поэтому обычной Gauss-Legendre квадратуры достаточно.
 */
static int rcs_integrate_triangle_phase(
    const double p0[3],
    const double p1[3],
    const double p2[3],
    double wave_number,
    const double observation_direction[3],
    int order,
    Complex *integral
)
{
    if (p0 == NULL ||
        p1 == NULL ||
        p2 == NULL ||
        observation_direction == NULL ||
        integral == NULL)
    {
        return -1;
    }

    if (order <= 0 || order > 16)
    {
        return -1;
    }


    double x[16];
    double w[16];


    if (gauss_legendre_rule(
            order,
            x,
            w
        ) != 0)
    {
        return -1;
    }


    integral->re =
        0.0;

    integral->im =
        0.0;


    for (int iu = 0; iu < order; ++iu)
    {
        double u =
            0.5 * (x[iu] + 1.0);

        double wu =
            0.5 * w[iu];


        for (int iv = 0; iv < order; ++iv)
        {
            double v =
                0.5 * (x[iv] + 1.0);

            double wv =
                0.5 * w[iv];


            double r[3];
            double jacobian = 0.0;


            rcs_triangle_map(
                p0,
                p1,
                p2,
                u,
                v,
                r,
                &jacobian
            );


            /*
             * Для exp(j*w*t) и исходящей волны
             * exp(-j*k*r) дальняя зона содержит:
             *
             *     exp(
             *         +j*k*r_hat.r'
             *     ).
             */
            double phase =
                wave_number
                * rcs_dot3(
                    observation_direction,
                    r
                );


            double weight =
                wu
                * wv
                * jacobian;


            integral->re +=
                weight
                * cos(phase);

            integral->im +=
                weight
                * sin(phase);
        }
    }


    return 0;
}


/*
 * ============================================================
 * ИНТЕГРАЛ ПО QUADPATCH
 * ============================================================
 *
 * Используем то же разбиение:
 *
 *     T0 = (0,1,2)
 *     T1 = (0,2,3).
 */
static int rcs_integrate_patch_phase(
    const QuadPatch *patch,
    double wave_number,
    const double observation_direction[3],
    int order,
    Complex *integral
)
{
    if (patch == NULL ||
        observation_direction == NULL ||
        integral == NULL)
    {
        return -1;
    }


    double p[4][3];


    for (int i = 0; i < 4; ++i)
    {
        rcs_patch_get_point(
            patch,
            i,
            p[i]
        );
    }


    Complex value0;
    Complex value1;


    if (rcs_integrate_triangle_phase(
            p[0],
            p[1],
            p[2],
            wave_number,
            observation_direction,
            order,
            &value0
        ) != 0)
    {
        return -1;
    }


    if (rcs_integrate_triangle_phase(
            p[0],
            p[2],
            p[3],
            wave_number,
            observation_direction,
            order,
            &value1
        ) != 0)
    {
        return -1;
    }


    integral->re =
        value0.re
        + value1.re;

    integral->im =
        value0.im
        + value1.im;


    return 0;
}


/*
 * ============================================================
 * ИНТЕГРАЛ ПО ВСЕЙ Pi_e^e
 * ============================================================
 *
 * Дуальная область ребра состоит из PatchList:
 *
 *     dual->edge_patches[e].
 *
 *
 * Поэтому:
 *
 *     G_e =
 *
 *     sum_p
 *
 *     integral_{patch p}
 *
 *     exp(j*k*r_hat.r') dS'.
 */
static int rcs_integrate_edge_region_phase(
    const DualMesh *dual,
    size_t edge,
    double wave_number,
    const double observation_direction[3],
    int order,
    Complex *integral
)
{
    if (dual == NULL ||
        observation_direction == NULL ||
        integral == NULL ||
        edge >= dual->n_edges)
    {
        return -1;
    }


    const PatchList *list =
        &dual->edge_patches[edge];


    integral->re =
        0.0;

    integral->im =
        0.0;


    for (size_t p = 0; p < list->count; ++p)
    {
        Complex patch_integral;


        if (rcs_integrate_patch_phase(
                &list->patches[p],
                wave_number,
                observation_direction,
                order,
                &patch_integral
            ) != 0)
        {
            return -1;
        }


        integral->re +=
            patch_integral.re;

        integral->im +=
            patch_integral.im;
    }


    return 0;
}


/*
 * ============================================================
 * ПОПЕРЕЧНЫЙ ВЕКТОР
 * ============================================================
 *
 *     T =
 *
 *     r_hat x (
 *         r_hat x F
 *     )
 *
 *
 * Используем:
 *
 *     a x (a x F)
 *
 *     = a(a.F) - F,
 *
 * при |a|=1.
 */
static void rcs_compute_transverse_vector(
    const double direction[3],
    const Complex F[3],
    Complex T[3]
)
{
    Complex dot;


    dot.re =
        direction[0] * F[0].re
        + direction[1] * F[1].re
        + direction[2] * F[2].re;


    dot.im =
        direction[0] * F[0].im
        + direction[1] * F[1].im
        + direction[2] * F[2].im;


    for (int c = 0; c < 3; ++c)
    {
        T[c].re =
            direction[c] * dot.re
            - F[c].re;

        T[c].im =
            direction[c] * dot.im
            - F[c].im;
    }
}


/*
 * ============================================================
 * |COMPLEX VECTOR|^2
 * ============================================================
 */
static double rcs_complex_vector_norm_squared(
    const Complex vector[3]
)
{
    double value =
        0.0;


    for (int c = 0; c < 3; ++c)
    {
        value +=
            vector[c].re * vector[c].re
            + vector[c].im * vector[c].im;
    }


    return value;
}


/*
 * ============================================================
 * OPTIONS
 * ============================================================
 */
void rcs_options_default(
    RcsOptions *options,
    int parallel_threads
)
{
    if (options == NULL)
    {
        return;
    }


    /*
     * Новый вариант теперь делаем основным.
     */
    options->integration =
        RCS_CURRENT_DUAL_PATCH_QUADRATURE;


    /*
     * Интеграл гладкий.
     *
     * Для ka=0.1 даже order=3-4 должен быть достаточен,
     * но 6 дает хороший запас и всё равно очень дешев.
     */
    options->surface_order =
        6;


    options->parallel_threads =
        parallel_threads;
}


/*
 * ============================================================
 * НАПРАВЛЕНИЕ НАБЛЮДЕНИЯ
 * ============================================================
 */
int rcs_build_observation_direction(
    double theta_deg,
    double phi_deg,
    double direction[3]
)
{
    if (direction == NULL ||
        !isfinite(theta_deg) ||
        !isfinite(phi_deg))
    {
        return -1;
    }


    double theta =
        theta_deg
        * M_PI
        / 180.0;


    double phi =
        phi_deg
        * M_PI
        / 180.0;


    direction[0] =
        sin(theta) * cos(phi);

    direction[1] =
        sin(theta) * sin(phi);

    direction[2] =
        cos(theta);


    double norm =
        rcs_norm3(
            direction
        );


    if (!isfinite(norm) ||
        norm <= 0.0)
    {
        return -1;
    }


    direction[0] /= norm;
    direction[1] /= norm;
    direction[2] /= norm;


    return 0;
}


/*
 * ============================================================
 * RADIATION VECTOR
 * ============================================================
 */
int rcs_compute_radiation_vector(
    const Mesh *mesh,
    const DualMesh *dual,
    const Complex *edge_current,
    double wave_number,
    const double observation_direction[3],
    const RcsOptions *options,
    Complex radiation_vector[3]
)
{
    if (mesh == NULL ||
        dual == NULL ||
        edge_current == NULL ||
        observation_direction == NULL ||
        options == NULL ||
        radiation_vector == NULL)
    {
        return -1;
    }


    if (mesh->n_edges == 0 ||
        mesh->edge_vectors == NULL)
    {
        fprintf(
            stderr,
            "ERROR: RCS requires edge vectors.\n"
        );

        return -1;
    }


    if (dual->n_edges != mesh->n_edges ||
        dual->edge_patches == NULL ||
        dual->edge_region_areas == NULL)
    {
        fprintf(
            stderr,
            "ERROR: RCS dual mesh does not match primal mesh.\n"
        );

        return -1;
    }


    if (!isfinite(wave_number) ||
        wave_number <= 0.0)
    {
        return -1;
    }


    if (options->parallel_threads < 1 ||
        (options->integration != RCS_CURRENT_MIDPOINT &&
         options->integration != RCS_CURRENT_DUAL_PATCH_QUADRATURE) ||
        !gauss_legendre_order_supported(options->surface_order)) {
        fprintf(stderr, "ERROR: invalid RCS integration or quadrature order.\n");
        return -1;
    }
    int failed = 0;
    double Fx_re = 0.0;
    double Fx_im = 0.0;

    double Fy_re = 0.0;
    double Fy_im = 0.0;

    double Fz_re = 0.0;
    double Fz_im = 0.0;


    #pragma omp parallel for \
        reduction(+:Fx_re,Fx_im,Fy_re,Fy_im,Fz_re,Fz_im) reduction(|:failed) \
        schedule(static) num_threads(options->parallel_threads)
    for (long long e_signed = 0;
         e_signed < (long long)mesh->n_edges;
         ++e_signed)
    {
        size_t e =
            (size_t)e_signed;


        const double *edge_vector =
            &mesh->edge_vectors[3 * e];


        double length =
            rcs_norm3(
                edge_vector
            );


        if (length <= 0.0)
        {
            continue;
        }


        /*
         * Единичный вектор ориентации ребра:
         *
         *     e_hat = l_vec / l.
         */
        double edge_direction[3] = {
            edge_vector[0] / length,
            edge_vector[1] / length,
            edge_vector[2] / length
        };


        Complex scalar;


        /*
         * ====================================================
         * MIDPOINT
         * ====================================================
         */
        if (options->integration ==
            RCS_CURRENT_MIDPOINT)
        {
            const double *center =
                &mesh->edge_centers[3 * e];


            double phase =
                wave_number
                * rcs_dot3(
                    observation_direction,
                    center
                );


            Complex exponential =
                rcs_complex_exp_j(
                    phase
                );


            /*
             * Старый предел:
             *
             *     F_e =
             *
             *     I_e * l_e * exp(j*phase).
             */
            scalar =
                rcs_complex_mul(
                    edge_current[e],
                    exponential
                );


            scalar.re *=
                length;

            scalar.im *=
                length;
        }


        /*
         * ====================================================
         * DUAL PATCH QUADRATURE
         * ====================================================
         */
        else
        {
            double area =
                dual->edge_region_areas[e];


            if (area <= 0.0)
            {
                continue;
            }


            /*
             * Эффективная ширина:
             *
             *     w_e =
             *
             *     S_e^e / l_e.
             *
             * [m]
             */
            double width =
                area / length;


            Complex G;


            if (rcs_integrate_edge_region_phase(
                    dual,
                    e,
                    wave_number,
                    observation_direction,
                    options->surface_order,
                    &G
                ) != 0)
            {
                failed = 1;
                continue;
            }


            /*
             * Поверхностная плотность тока:
             *
             *     J_e =
             *
             *     I_e / w_e
             *     * e_hat.
             *
             *
             * Поэтому скалярная часть:
             *
             *     I_e / w_e * G_e.
             */
            Complex current_density_factor =
                edge_current[e];


            current_density_factor.re /=
                width;

            current_density_factor.im /=
                width;


            scalar =
                rcs_complex_mul(
                    current_density_factor,
                    G
                );
        }


        /*
         * Добавляем:
         *
         *     scalar * e_hat
         *
         * к F.
         */
        Fx_re +=
            scalar.re
            * edge_direction[0];

        Fx_im +=
            scalar.im
            * edge_direction[0];


        Fy_re +=
            scalar.re
            * edge_direction[1];

        Fy_im +=
            scalar.im
            * edge_direction[1];


        Fz_re +=
            scalar.re
            * edge_direction[2];

        Fz_im +=
            scalar.im
            * edge_direction[2];
    }


    if (failed) {
        fprintf(stderr, "ERROR: RCS surface quadrature failed.\n");
        return -1;
    }

    radiation_vector[0].re = Fx_re;
    radiation_vector[0].im = Fx_im;

    radiation_vector[1].re = Fy_re;
    radiation_vector[1].im = Fy_im;

    radiation_vector[2].re = Fz_re;
    radiation_vector[2].im = Fz_im;


    return 0;
}


/*
 * ============================================================
 * RCS В ОДНОМ НАПРАВЛЕНИИ
 * ============================================================
 */
int rcs_compute_direction(
    const Mesh *mesh,
    const DualMesh *dual,
    const IncidentField *incident_field,
    const HarmonicSolution *solution,
    const RcsOptions *options,
    double theta_deg,
    double phi_deg,
    RcsResult *result
)
{
    if (mesh == NULL ||
        dual == NULL ||
        incident_field == NULL ||
        solution == NULL ||
        options == NULL ||
        result == NULL)
    {
        return -1;
    }


    if (solution->edge_current == NULL ||
        solution->n_edges != mesh->n_edges)
    {
        fprintf(
            stderr,
            "ERROR: RCS solution does not match mesh.\n"
        );

        return -1;
    }


    if (!isfinite(incident_field->E_amp) ||
        incident_field->E_amp <= 0.0)
    {
        fprintf(
            stderr,
            "ERROR: RCS requires positive incident field amplitude.\n"
        );

        return -1;
    }


    if (!isfinite(incident_field->wave_number) ||
        incident_field->wave_number <= 0.0)
    {
        fprintf(
            stderr,
            "ERROR: RCS requires positive wave number.\n"
        );

        return -1;
    }


    memset(
        result,
        0,
        sizeof(*result)
    );


    result->theta_deg =
        theta_deg;

    result->phi_deg =
        phi_deg;


    if (rcs_build_observation_direction(
            theta_deg,
            phi_deg,
            result->direction
        ) != 0)
    {
        return -1;
    }


    if (rcs_compute_radiation_vector(
            mesh,
            dual,
            solution->edge_current,
            incident_field->wave_number,
            result->direction,
            options,
            result->radiation_vector
        ) != 0)
    {
        return -1;
    }


    rcs_compute_transverse_vector(
        result->direction,
        result->radiation_vector,
        result->transverse_vector
    );


    double transverse_norm_squared =
        rcs_complex_vector_norm_squared(
            result->transverse_vector
        );


    double eta0 =
        sqrt(
            RCS_MU0
            / RCS_EPS0
        );


    double k =
        incident_field->wave_number;


    double E0 =
        incident_field->E_amp;


    /*
     *     sigma =
     *
     *     k^2 eta0^2
     *     ------------
     *       4 pi E0^2
     *
     *     |T|^2.
     */
    result->sigma =
        k * k
        * eta0 * eta0
        * transverse_norm_squared
        / (
            4.0
            * M_PI
            * E0 * E0
        );


    if (result->sigma > DBL_MIN)
    {
        result->sigma_dbsm =
            10.0
            * log10(
                result->sigma
            );
    }
    else
    {
        result->sigma_dbsm =
            -INFINITY;
    }


    return 0;
}


/*
 * ============================================================
 * ДАЛЬНЕЕ E_sca
 * ============================================================
 */
int rcs_compute_far_electric_field(
    const RcsResult *result,
    double wave_number,
    double distance,
    Complex E_sca[3]
)
{
    if (result == NULL ||
        E_sca == NULL ||
        !isfinite(wave_number) ||
        !isfinite(distance) ||
        wave_number <= 0.0 ||
        distance <= 0.0)
    {
        return -1;
    }


    double eta0 =
        sqrt(
            RCS_MU0
            / RCS_EPS0
        );


    double kr =
        wave_number
        * distance;


    double amplitude =
        wave_number
        * eta0
        / (
            4.0
            * M_PI
            * distance
        );


    /*
     *     j exp(-jkr)
     *
     *     = sin(kr) + j cos(kr).
     */
    Complex coefficient;

    coefficient.re =
        amplitude
        * sin(kr);

    coefficient.im =
        amplitude
        * cos(kr);


    for (int c = 0; c < 3; ++c)
    {
        E_sca[c] =
            rcs_complex_mul(
                coefficient,
                result->transverse_vector[c]
            );
    }


    return 0;
}


/*
 * ============================================================
 * ВЫВОД
 * ============================================================
 */
void rcs_print_result(
    const RcsResult *result
)
{
    if (result == NULL)
    {
        return;
    }


    printf("\n");
    printf("========================================\n");
    printf("RADAR CROSS SECTION\n");
    printf("========================================\n");

    printf(
        "Observation theta : %.6f deg\n",
        result->theta_deg
    );

    printf(
        "Observation phi   : %.6f deg\n",
        result->phi_deg
    );

    printf(
        "Direction         : %.9e %.9e %.9e\n",
        result->direction[0],
        result->direction[1],
        result->direction[2]
    );


    printf("\n");

    printf(
        "F_x               : %.9e %+.9ej A*m\n",
        result->radiation_vector[0].re,
        result->radiation_vector[0].im
    );

    printf(
        "F_y               : %.9e %+.9ej A*m\n",
        result->radiation_vector[1].re,
        result->radiation_vector[1].im
    );

    printf(
        "F_z               : %.9e %+.9ej A*m\n",
        result->radiation_vector[2].re,
        result->radiation_vector[2].im
    );


    printf("\n");

    printf(
        "RCS sigma         : %.9e m^2\n",
        result->sigma
    );

    printf(
        "RCS               : %.9f dBsm\n",
        result->sigma_dbsm
    );

    printf("========================================\n");
}