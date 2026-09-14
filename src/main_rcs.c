#include "mesh.h"
#include "topology.h"
#include "dual_mesh.h"

#include "matrix.h"

#include "incidence.h"
#include "inductance.h"
#include "potential.h"
#include "resistance.h"

#include "harmonic.h"

#include "incident_field.h"
#include "scattering.h"
#include "rcs.h"

#include <errno.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


/*
 * ============================================================
 * ФИЗИЧЕСКИЕ КОНСТАНТЫ
 * ============================================================
 */

static const double C0 = 299792458.0;


/*
 * ============================================================
 * СПОСОБ ЗАДАНИЯ ЧАСТОТЫ
 * ============================================================
 */
typedef enum
{
    RCS_FREQUENCY_FROM_KA = 0,
    RCS_FREQUENCY_DIRECT,
    RCS_FREQUENCY_FROM_K

} RcsFrequencyMode;


/*
 * ============================================================
 * ПАРАМЕТРЫ ЗАДАЧИ
 * ============================================================
 */
typedef struct
{
    /*
     * Сетка.
     */
    const char *mesh_file;

    /*
     * OpenMP.
     */
    int parallel_threads;


    /*
     * --------------------------------------------------------
     * ЭЛЕКТРИЧЕСКИЙ РАЗМЕР
     * --------------------------------------------------------
     *
     * Можно задать один из вариантов:
     *
     *     --frequency f
     *
     * или:
     *
     *     --k k
     *
     * или:
     *
     *     --ka ka --a a.
     */
    RcsFrequencyMode frequency_mode;

    double frequency;
    double wave_number;

    double ka;
    double characteristic_size;


    /*
     * --------------------------------------------------------
     * ПАДАЮЩАЯ ВОЛНА
     * --------------------------------------------------------
     */

    double E_amplitude;

    double incident_theta;
    double incident_phi;

    double incident_phase;

    IncidentPolarization polarization;


    /*
     * --------------------------------------------------------
     * RCS SWEEP
     * --------------------------------------------------------
     */

    double observation_theta;

    double phi_start;
    double phi_end;
    double phi_step;


    /*
     * --------------------------------------------------------
     * ИНТЕГРИРОВАНИЕ ДАЛЬНЕЙ ЗОНЫ
     * --------------------------------------------------------
     */

    RcsIntegration rcs_integration;

    int rcs_surface_order;


    /*
     * --------------------------------------------------------
     * ВЫХОДНОЙ ФАЙЛ
     * --------------------------------------------------------
     */

    const char *output_file;

} RcsRunConfig;


/*
 * ============================================================
 * DEFAULT
 * ============================================================
 *
 * По умолчанию:
 *
 *     ka = 0.1
 *     a  = 1 m.
 *
 *
 * Поэтому:
 *
 *     k = ka/a = 0.1 rad/m
 *
 * и:
 *
 *     f = k*c0/(2*pi).
 */
static void rcs_run_config_default(
    RcsRunConfig *config
)
{
    if (config == NULL)
    {
        return;
    }

    config->mesh_file = "meshes/body_closed.msh";

    config->parallel_threads = 12;


    /*
     * Частота по умолчанию определяется через ka.
     */
    config->frequency_mode = RCS_FREQUENCY_FROM_KA;

    config->frequency = 0.0;
    config->wave_number = 0.0;

    config->ka = 0.1;
    config->characteristic_size = 1.0;


    /*
     * Падающая волна:
     *
     *     theta = 90 deg
     *     phi   = 0 deg
     *
     * при нашей конвенции:
     *
     *     k_hat_inc = (-1,0,0).
     */
    config->E_amplitude = 1.0;

    config->incident_theta = 90.0;
    config->incident_phi = 0.0;

    config->incident_phase = 0.0;

    config->polarization =
        INCIDENT_POLARIZATION_HORIZONTAL;


    /*
     * RCS в плоскости xy.
     */
    config->observation_theta = 90.0;

    config->phi_start = 0.0;
    config->phi_end = 360.0;
    config->phi_step = 1.0;


    /*
     * Полное интегрирование по Pi_e^e.
     */
    config->rcs_integration =
        RCS_CURRENT_DUAL_PATCH_QUADRATURE;

    config->rcs_surface_order = 6;


    config->output_file =
        "results/rcs.csv";
}


/*
 * ============================================================
 * HELP
 * ============================================================
 */
static void print_help(
    const char *program
)
{
    printf("\n");
    printf("PEEC RCS calculation\n");
    printf("\n");

    printf("Usage:\n");
    printf("\n");
    printf("  %s [options]\n", program);
    printf("\n");

    printf("Mesh and parallelism:\n");
    printf("  --mesh FILE              Gmsh mesh file\n");
    printf("  --parallel N             OpenMP threads\n");
    printf("\n");

    printf("Frequency specification:\n");
    printf("  --frequency F            Frequency [Hz]\n");
    printf("  --k K                    Wave number [rad/m]\n");
    printf("  --ka KA                  Electrical size k*a\n");
    printf("  --a A                    Characteristic body size [m]\n");
    printf("\n");

    printf("Incident wave:\n");
    printf("  --E VALUE                Incident E amplitude [V/m]\n");
    printf("  --theta VALUE            Incident theta [deg]\n");
    printf("  --phi VALUE              Incident phi [deg]\n");
    printf("  --phase VALUE            Initial phase [rad]\n");
    printf("  --polarization horizontal|vertical\n");
    printf("\n");

    printf("RCS observation sweep:\n");
    printf("  --obs-theta VALUE        Observation theta [deg]\n");
    printf("  --phi-start VALUE        Sweep start phi [deg]\n");
    printf("  --phi-end VALUE          Sweep end phi [deg]\n");
    printf("  --phi-step VALUE         Sweep step phi [deg]\n");
    printf("\n");

    printf("RCS integration:\n");
    printf("  --rcs-midpoint           Old midpoint approximation\n");
    printf("  --rcs-dual               Integrate over dual edge regions\n");
    printf("  --rcs-order N            Gauss-Legendre order\n");
    printf("\n");

    printf("Output:\n");
    printf("  --output FILE            CSV file\n");
    printf("\n");

    printf("Other:\n");
    printf("  --help                   Show this help\n");
    printf("\n");

    printf("Examples:\n");
    printf("\n");

    printf(
        "  %s --ka 0.1 --a 1.0\n",
        program
    );

    printf(
        "  %s --frequency 50e6\n",
        program
    );

    printf(
        "  %s --k 0.5 --a 1.0\n",
        program
    );

    printf(
        "  %s --ka 1 --a 1 --theta 90 --phi 45 "
        "--polarization horizontal\n",
        program
    );

    printf("\n");
}


/*
 * ============================================================
 * ЧТЕНИЕ DOUBLE
 * ============================================================
 */
static int parse_double_value(
    const char *text,
    double *value
)
{
    if (text == NULL || value == NULL)
    {
        return -1;
    }

    char *end = NULL;

    errno = 0;

    double result =
        strtod(
            text,
            &end
        );

    if (errno != 0 ||
        end == text ||
        *end != '\0' ||
        !isfinite(result))
    {
        return -1;
    }

    *value =
        result;

    return 0;
}


/*
 * ============================================================
 * ЧТЕНИЕ INTEGER
 * ============================================================
 */
static int parse_int_value(
    const char *text,
    int *value
)
{
    if (text == NULL || value == NULL)
    {
        return -1;
    }

    char *end = NULL;

    errno = 0;

    long result =
        strtol(
            text,
            &end,
            10
        );

    if (errno != 0 ||
        end == text ||
        *end != '\0' ||
        result < 1 ||
        result > 1000000)
    {
        return -1;
    }

    *value =
        (int)result;

    return 0;
}


/*
 * ============================================================
 * COMMAND LINE
 * ============================================================
 */
static int rcs_parse_arguments(
    RcsRunConfig *config,
    int argc,
    char **argv
)
{
    if (config == NULL)
    {
        return -1;
    }


    /*
     * Считаем, сколько разных способов задания частоты
     * пользователь явно указал.
     */
    int frequency_given = 0;
    int k_given = 0;
    int ka_given = 0;


    for (int i = 1; i < argc; ++i)
    {
        /*
         * ----------------------------------------------------
         * HELP
         * ----------------------------------------------------
         */
        if (strcmp(argv[i], "--help") == 0 ||
            strcmp(argv[i], "-h") == 0)
        {
            print_help(
                argv[0]
            );

            return 1;
        }


        /*
         * ----------------------------------------------------
         * MESH
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--mesh") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --mesh requires a file.\n"
                );

                return -1;
            }

            config->mesh_file =
                argv[++i];
        }


        /*
         * ----------------------------------------------------
         * OPENMP
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--parallel") == 0)
        {
            if (i + 1 >= argc ||
                parse_int_value(
                    argv[++i],
                    &config->parallel_threads
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --parallel value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * FREQUENCY [Hz]
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--frequency") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->frequency
                ) != 0 ||
                config->frequency <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: --frequency must be positive [Hz].\n"
                );

                return -1;
            }

            ++frequency_given;
        }


        /*
         * ----------------------------------------------------
         * WAVE NUMBER k [rad/m]
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--k") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->wave_number
                ) != 0 ||
                config->wave_number <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: --k must be positive [rad/m].\n"
                );

                return -1;
            }

            ++k_given;
        }


        /*
         * ----------------------------------------------------
         * ELECTRICAL SIZE ka
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--ka") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->ka
                ) != 0 ||
                config->ka <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: --ka must be positive.\n"
                );

                return -1;
            }

            ++ka_given;
        }


        /*
         * ----------------------------------------------------
         * CHARACTERISTIC SIZE a [m]
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--a") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->characteristic_size
                ) != 0 ||
                config->characteristic_size <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: --a must be positive [m].\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * INCIDENT FIELD AMPLITUDE
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--E") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->E_amplitude
                ) != 0 ||
                config->E_amplitude <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: --E must be positive [V/m].\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * INCIDENT THETA
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--theta") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->incident_theta
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --theta value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * INCIDENT PHI
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--phi") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->incident_phi
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phi value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * INITIAL PHASE
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--phase") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->incident_phase
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phase value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * POLARIZATION
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--polarization") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --polarization requires horizontal or vertical.\n"
                );

                return -1;
            }

            const char *value =
                argv[++i];

            if (strcmp(value, "horizontal") == 0 ||
                strcmp(value, "h") == 0)
            {
                config->polarization =
                    INCIDENT_POLARIZATION_HORIZONTAL;
            }
            else if (strcmp(value, "vertical") == 0 ||
                     strcmp(value, "v") == 0)
            {
                config->polarization =
                    INCIDENT_POLARIZATION_VERTICAL;
            }
            else
            {
                fprintf(
                    stderr,
                    "ERROR: unknown polarization '%s'.\n",
                    value
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * OBSERVATION THETA
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--obs-theta") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->observation_theta
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --obs-theta value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * PHI START
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--phi-start") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->phi_start
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phi-start value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * PHI END
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--phi-end") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->phi_end
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phi-end value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * PHI STEP
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--phi-step") == 0)
        {
            if (i + 1 >= argc ||
                parse_double_value(
                    argv[++i],
                    &config->phi_step
                ) != 0 ||
                config->phi_step <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: --phi-step must be positive.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * OLD RCS MIDPOINT
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--rcs-midpoint") == 0)
        {
            config->rcs_integration =
                RCS_CURRENT_MIDPOINT;
        }


        /*
         * ----------------------------------------------------
         * DUAL PATCH RCS
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--rcs-dual") == 0)
        {
            config->rcs_integration =
                RCS_CURRENT_DUAL_PATCH_QUADRATURE;
        }


        /*
         * ----------------------------------------------------
         * SURFACE QUADRATURE ORDER
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--rcs-order") == 0)
        {
            if (i + 1 >= argc ||
                parse_int_value(
                    argv[++i],
                    &config->rcs_surface_order
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --rcs-order value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * OUTPUT
         * ----------------------------------------------------
         */
        else if (strcmp(argv[i], "--output") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --output requires a filename.\n"
                );

                return -1;
            }

            config->output_file =
                argv[++i];
        }


        /*
         * ----------------------------------------------------
         * UNKNOWN OPTION
         * ----------------------------------------------------
         */
        else
        {
            fprintf(
                stderr,
                "ERROR: unknown argument '%s'.\n",
                argv[i]
            );

            fprintf(
                stderr,
                "Use --help for available options.\n"
            );

            return -1;
        }
    }


    /*
     * Нельзя одновременно явно задавать:
     *
     *     frequency,
     *     k,
     *     ka.
     */
    int modes =
        (frequency_given > 0)
        + (k_given > 0)
        + (ka_given > 0);

    if (modes > 1)
    {
        fprintf(
            stderr,
            "ERROR: use only one of --frequency, --k or --ka.\n"
        );

        return -1;
    }


    if (frequency_given)
    {
        config->frequency_mode =
            RCS_FREQUENCY_DIRECT;
    }
    else if (k_given)
    {
        config->frequency_mode =
            RCS_FREQUENCY_FROM_K;
    }
    else
    {
        /*
         * Если ничего не задано или явно задан --ka,
         * используем ka.
         */
        config->frequency_mode =
            RCS_FREQUENCY_FROM_KA;
    }


    return 0;
}


/*
 * ============================================================
 * РАСЧЕТ f, k И ka
 * ============================================================
 */
static int rcs_resolve_frequency(
    RcsRunConfig *config
)
{
    if (config == NULL ||
        config->characteristic_size <= 0.0)
    {
        return -1;
    }


    /*
     * --------------------------------------------------------
     * f задан непосредственно
     * --------------------------------------------------------
     */
    if (config->frequency_mode ==
        RCS_FREQUENCY_DIRECT)
    {
        config->wave_number =
            2.0
            * M_PI
            * config->frequency
            / C0;

        config->ka =
            config->wave_number
            * config->characteristic_size;
    }


    /*
     * --------------------------------------------------------
     * k задан непосредственно
     * --------------------------------------------------------
     */
    else if (config->frequency_mode ==
             RCS_FREQUENCY_FROM_K)
    {
        config->frequency =
            config->wave_number
            * C0
            / (
                2.0
                * M_PI
            );

        config->ka =
            config->wave_number
            * config->characteristic_size;
    }


    /*
     * --------------------------------------------------------
     * ka задано
     * --------------------------------------------------------
     *
     *     ka = k*a
     *
     * поэтому:
     *
     *     k = ka/a.
     *
     *
     * Затем:
     *
     *     f = k*c0/(2*pi).
     */
    else
    {
        config->wave_number =
            config->ka
            / config->characteristic_size;

        config->frequency =
            config->wave_number
            * C0
            / (
                2.0
                * M_PI
            );
    }


    if (!isfinite(config->frequency) ||
        !isfinite(config->wave_number) ||
        !isfinite(config->ka) ||
        config->frequency <= 0.0 ||
        config->wave_number <= 0.0 ||
        config->ka <= 0.0)
    {
        return -1;
    }


    return 0;
}


/*
 * ============================================================
 * ПЕЧАТЬ ПАРАМЕТРОВ
 * ============================================================
 */
static void rcs_print_run_config(
    const RcsRunConfig *config
)
{
    if (config == NULL)
    {
        return;
    }


    double wavelength =
        2.0
        * M_PI
        / config->wave_number;


    printf("\n");
    printf("============================================================\n");
    printf("PEEC RCS CALCULATION\n");
    printf("============================================================\n");

    printf(
        "Mesh               : %s\n",
        config->mesh_file
    );

    printf(
        "Parallel           : %d OpenMP thread(s)\n",
        config->parallel_threads
    );

    printf(
        "Physics            : QUASISTATIC PEEC\n"
    );


    printf("\n");
    printf("Electrical size\n");
    printf("------------------------------------------------------------\n");

    printf(
        "a                  : %.9e m\n",
        config->characteristic_size
    );

    printf(
        "ka                 : %.9e\n",
        config->ka
    );

    printf(
        "k                  : %.9e rad/m\n",
        config->wave_number
    );

    printf(
        "Frequency          : %.9e Hz\n",
        config->frequency
    );

    printf(
        "Wavelength         : %.9e m\n",
        wavelength
    );


    printf("\n");
    printf("Incident wave\n");
    printf("------------------------------------------------------------\n");

    printf(
        "E amplitude        : %.9e V/m\n",
        config->E_amplitude
    );

    printf(
        "Theta              : %.6f deg\n",
        config->incident_theta
    );

    printf(
        "Phi                : %.6f deg\n",
        config->incident_phi
    );

    printf(
        "Phase              : %.9e rad\n",
        config->incident_phase
    );

    printf(
        "Polarization       : %s\n",
        config->polarization ==
        INCIDENT_POLARIZATION_HORIZONTAL
            ? "HORIZONTAL"
            : "VERTICAL"
    );


    printf("\n");
    printf("RCS sweep\n");
    printf("------------------------------------------------------------\n");

    printf(
        "Observation theta  : %.6f deg\n",
        config->observation_theta
    );

    printf(
        "Phi                : %.6f ... %.6f deg\n",
        config->phi_start,
        config->phi_end
    );

    printf(
        "Phi step           : %.6f deg\n",
        config->phi_step
    );

    printf(
        "Integration        : %s\n",
        config->rcs_integration ==
        RCS_CURRENT_DUAL_PATCH_QUADRATURE
            ? "DUAL PATCH QUADRATURE"
            : "MIDPOINT"
    );

    if (config->rcs_integration ==
        RCS_CURRENT_DUAL_PATCH_QUADRATURE)
    {
        printf(
            "Surface order      : %d\n",
            config->rcs_surface_order
        );
    }

    printf(
        "Output             : %s\n",
        config->output_file
    );

    printf("============================================================\n");
}


/*
 * ============================================================
 * MAX |COMPLEX|
 * ============================================================
 */
static double complex_array_max_abs(
    const Complex *values,
    size_t count
)
{
    if (values == NULL || count == 0)
    {
        return 0.0;
    }


    double maximum =
        0.0;


    for (size_t i = 0; i < count; ++i)
    {
        double magnitude =
            hypot(
                values[i].re,
                values[i].im
            );

        if (magnitude > maximum)
        {
            maximum =
                magnitude;
        }
    }


    return maximum;
}


/*
 * ============================================================
 * results/
 * ============================================================
 */
static int ensure_results_directory(void)
{
    if (mkdir("results", 0755) == 0)
    {
        return 0;
    }

    if (errno == EEXIST)
    {
        return 0;
    }

    fprintf(
        stderr,
        "ERROR: cannot create results directory: %s\n",
        strerror(errno)
    );

    return -1;
}


/*
 * ============================================================
 * RCS SWEEP -> CSV
 * ============================================================
 */
static int save_rcs_sweep(
    const RcsRunConfig *config,
    const Mesh *mesh,
    const DualMesh *dual,
    const IncidentField *field,
    const HarmonicSolution *solution,
    const RcsOptions *options
)
{
    if (config == NULL ||
        mesh == NULL ||
        dual == NULL ||
        field == NULL ||
        solution == NULL ||
        options == NULL)
    {
        return -1;
    }


    FILE *file =
        fopen(
            config->output_file,
            "w"
        );


    if (file == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot open '%s'.\n",
            config->output_file
        );

        return -1;
    }


    /*
     * Дополнительно сохраняем параметры задачи прямо
     * в CSV как обычные дополнительные колонки.
     *
     * Это удобно: файл потом однозначно знает, при каком
     * ka, k, f и a был получен.
     */
    fprintf(
        file,
        "theta_deg,"
        "phi_deg,"
        "sigma_m2,"
        "sigma_dbsm,"
        "frequency_hz,"
        "k_rad_m,"
        "ka,"
        "a_m\n"
    );


    size_t point_count =
        0;


    double maximum_sigma =
        -1.0;

    double maximum_dbsm =
        -INFINITY;

    double maximum_phi =
        0.0;


    /*
     * Сохраняем отдельно значение в направлении
     * обратного рассеяния.
     *
     * Оно будет вычислено ниже по геометрии падающей волны.
     */
    double backscatter_sigma =
        NAN;

    double backscatter_dbsm =
        NAN;


    /*
     * Направление обратно к источнику:
     *
     *     r_hat_back = -k_hat_inc.
     *
     *
     * Переводим его в сферический phi в плоскости sweep,
     * если observation theta совпадает с его theta.
     */
    double back_x =
        -field->k_hat[0];

    double back_y =
        -field->k_hat[1];

    double back_z =
        -field->k_hat[2];


    double back_theta =
        acos(
            fmax(
                -1.0,
                fmin(
                    1.0,
                    back_z
                )
            )
        )
        * 180.0
        / M_PI;


    double back_phi =
        atan2(
            back_y,
            back_x
        )
        * 180.0
        / M_PI;


    if (back_phi < 0.0)
    {
        back_phi +=
            360.0;
    }


    /*
     * Сначала считаем точный backscatter независимо
     * от sweep.
     */
    RcsResult backscatter;


    if (rcs_compute_direction(
            mesh,
            dual,
            field,
            solution,
            options,
            back_theta,
            back_phi,
            &backscatter
        ) != 0)
    {
        fclose(file);

        fprintf(
            stderr,
            "ERROR: backscatter RCS calculation failed.\n"
        );

        return -1;
    }


    backscatter_sigma =
        backscatter.sigma;

    backscatter_dbsm =
        backscatter.sigma_dbsm;


    /*
     * Основной angular sweep.
     */
    for (double phi = config->phi_start;
         phi <= config->phi_end + 0.5 * config->phi_step;
         phi += config->phi_step)
    {
        RcsResult result;


        if (rcs_compute_direction(
                mesh,
                dual,
                field,
                solution,
                options,
                config->observation_theta,
                phi,
                &result
            ) != 0)
        {
            fclose(file);

            fprintf(
                stderr,
                "ERROR: RCS calculation failed at phi = %.6f deg.\n",
                phi
            );

            return -1;
        }


        fprintf(
            file,
            "%.9f,"
            "%.9f,"
            "%.16e,"
            "%.16e,"
            "%.16e,"
            "%.16e,"
            "%.16e,"
            "%.16e\n",
            result.theta_deg,
            result.phi_deg,
            result.sigma,
            result.sigma_dbsm,
            config->frequency,
            config->wave_number,
            config->ka,
            config->characteristic_size
        );


        if (result.sigma > maximum_sigma)
        {
            maximum_sigma =
                result.sigma;

            maximum_dbsm =
                result.sigma_dbsm;

            maximum_phi =
                phi;
        }


        ++point_count;
    }


    fclose(file);


    printf("\n");
    printf("============================================================\n");
    printf("RCS RESULTS\n");
    printf("============================================================\n");

    printf(
        "Backscatter theta  : %.6f deg\n",
        back_theta
    );

    printf(
        "Backscatter phi    : %.6f deg\n",
        back_phi
    );

    printf(
        "Backscatter sigma  : %.9e m^2\n",
        backscatter_sigma
    );

    printf(
        "Backscatter RCS    : %.9f dBsm\n",
        backscatter_dbsm
    );


    printf("\n");

    printf(
        "Maximum sweep RCS  : %.9e m^2\n",
        maximum_sigma
    );

    printf(
        "Maximum sweep RCS  : %.9f dBsm\n",
        maximum_dbsm
    );

    printf(
        "Maximum at phi     : %.6f deg\n",
        maximum_phi
    );


    printf("\n");

    printf(
        "Sweep points       : %zu\n",
        point_count
    );

    printf(
        "Saved              : %s\n",
        config->output_file
    );

    printf("============================================================\n");


    return 0;
}


/*
 * ============================================================
 * MAIN
 * ============================================================
 */
int main(
    int argc,
    char **argv
)
{
    int exit_code =
        EXIT_FAILURE;


    /*
     * ========================================================
     * CONFIG
     * ========================================================
     */
    RcsRunConfig config;

    rcs_run_config_default(
        &config
    );


    int parse_status =
        rcs_parse_arguments(
            &config,
            argc,
            argv
        );


    if (parse_status > 0)
    {
        return EXIT_SUCCESS;
    }

    if (parse_status < 0)
    {
        return EXIT_FAILURE;
    }


    if (rcs_resolve_frequency(
            &config
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot determine frequency.\n"
        );

        return EXIT_FAILURE;
    }


    rcs_print_run_config(
        &config
    );


    /*
     * ========================================================
     * OPENMP
     * ========================================================
     */
    omp_set_dynamic(0);

    omp_set_num_threads(
        config.parallel_threads
    );


    /*
     * ========================================================
     * ОСНОВНЫЕ СТРУКТУРЫ
     * ========================================================
     */
    Mesh mesh;
    DualMesh dual;

    IncidenceMatrix incidence;

    InductanceMatrix L;
    PotentialMatrix P;

    Matrix PAT;

    ResistanceMatrix R;

    HarmonicSystem harmonic_system;
    HarmonicSolution solution;


    Complex *U_inc =
        NULL;


    mesh_init(
        &mesh
    );

    dual_mesh_init(
        &dual
    );

    incidence_init(
        &incidence
    );

    inductance_init(
        &L
    );

    potential_init(
        &P
    );

    matrix_init(
        &PAT
    );

    resistance_init(
        &R
    );

    harmonic_system_init(
        &harmonic_system
    );

    harmonic_solution_init(
        &solution
    );


    /*
     * ========================================================
     * 1. MESH
     * ========================================================
     */
    printf("\nReading Gmsh mesh ...\n");


    if (mesh_load_gmsh(
            config.mesh_file,
            &mesh
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot read mesh.\n"
        );

        goto cleanup;
    }


    printf("Building topology ...\n");


    if (mesh_build_topology(
            &mesh
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: topology construction failed.\n"
        );

        goto cleanup;
    }


    mesh_print_info(
        &mesh
    );


    /*
     * ========================================================
     * 2. DUAL MESH
     * ========================================================
     */
    printf("\nBuilding dual mesh ...\n");


    if (dual_mesh_build(
            &mesh,
            &dual
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: dual mesh construction failed.\n"
        );

        goto cleanup;
    }


    dual_mesh_print_info(
        &dual
    );


    /*
     * ========================================================
     * 3. INCIDENCE MATRIX A
     * ========================================================
     */
    printf("\nBuilding incidence matrix ...\n");


    if (incidence_build(
            &mesh,
            &incidence
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: incidence matrix construction failed.\n"
        );

        goto cleanup;
    }


    /*
     * ========================================================
     * 4. L
     * ========================================================
     */
    InductanceOptions inductance_options;


    inductance_options_default(
        &inductance_options,
        config.parallel_threads
    );


    printf("\nComputing inductance matrix L ...\n");


    double L_start =
        omp_get_wtime();


    if (inductance_compute_matrix(
            &mesh,
            &dual,
            &inductance_options,
            &L
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: inductance calculation failed.\n"
        );

        goto cleanup;
    }


    double L_time =
        omp_get_wtime()
        - L_start;


    printf(
        "L computation time : %.6f s\n",
        L_time
    );


    /*
     * ========================================================
     * 5. P
     * ========================================================
     */
    PotentialOptions potential_options;


    potential_options_default(
        &potential_options,
        config.parallel_threads
    );


    printf("\nComputing potential matrix P ...\n");


    double P_start =
        omp_get_wtime();


    if (potential_compute_matrix(
            &mesh,
            &dual,
            &potential_options,
            &P
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: potential calculation failed.\n"
        );

        goto cleanup;
    }


    double P_time =
        omp_get_wtime()
        - P_start;


    printf(
        "P computation time : %.6f s\n",
        P_time
    );


    /*
     * ========================================================
     * 6. P A^T
     * ========================================================
     */
    printf("\nBuilding P A^T ...\n");


    if (incidence_build_PAT(
            &mesh,
            &P,
            &PAT,
            config.parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: P A^T construction failed.\n"
        );

        goto cleanup;
    }


    /*
     * ========================================================
     * 7. R
     * ========================================================
     */
    ResistanceOptions resistance_options;


    resistance_options_default(
        &resistance_options,
        config.parallel_threads
    );


    resistance_options.model =
        PEEC_CONDUCTOR_PEC;


    printf("\nBuilding resistance matrix R ...\n");


    if (resistance_compute_matrix(
            &mesh,
            &dual,
            &resistance_options,
            &R
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: resistance calculation failed.\n"
        );

        goto cleanup;
    }


    /*
     * ========================================================
     * 8. INCIDENT FIELD
     * ========================================================
     */
    IncidentField field;


    incident_field_init(
        &field
    );


    if (incident_field_configure_harmonic(
            &field,
            INCIDENT_WAVE_COMPLEX,
            config.polarization,
            config.E_amplitude,
            config.frequency,
            config.incident_theta,
            config.incident_phi,
            config.incident_phase
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: incident field configuration failed.\n"
        );

        goto cleanup;
    }


    incident_field_print_info(
        &field
    );


    /*
     * Проверяем, что k, рассчитанный IncidentField,
     * совпадает с нашим k.
     */
    double k_error =
        fabs(
            field.wave_number
            - config.wave_number
        );


    double k_scale =
        fmax(
            1.0,
            fabs(config.wave_number)
        );


    if (k_error / k_scale > 1.0e-10)
    {
        fprintf(
            stderr,
            "ERROR: inconsistent wave number.\n"
        );

        goto cleanup;
    }


    /*
     * ========================================================
     * 9. U_inc
     * ========================================================
     */
    ScatteringOptions scattering_options;


    scattering_options_default(
        &scattering_options,
        config.parallel_threads
    );


    /*
     * Для гармонической плоской волны:
     *
     *     U_e =
     *
     *     integral_edge
     *
     *     E_inc . dl.
     */
    scattering_options.integration =
        SCATTERING_EDGE_ANALYTIC_PLANE_WAVE;


    U_inc =
        malloc(
            mesh.n_edges
            * sizeof(Complex)
        );


    if (U_inc == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate U_inc.\n"
        );

        goto cleanup;
    }


    if (scattering_build_harmonic_excitation(
            &mesh,
            &field,
            &scattering_options,
            U_inc
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: incident excitation calculation failed.\n"
        );

        goto cleanup;
    }


    scattering_print_harmonic_excitation_info(
        U_inc,
        mesh.n_edges
    );


    /*
     * ========================================================
     * 10. HARMONIC SYSTEM
     * ========================================================
     *
     * Полная классическая система:
     *
     *     [ R+jwL      A  ] [I]   [-U_inc]
     *     [                ] [ ] = [      ]
     *     [-PA^T      jwIv] [V]   [   0  ].
     */
    printf("\nBuilding harmonic PEEC system ...\n");


    double matrix_start =
        omp_get_wtime();


    if (harmonic_build_system(
            &L,
            &R,
            &incidence,
            &PAT,
            config.frequency,
            config.parallel_threads,
            &harmonic_system
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic matrix construction failed.\n"
        );

        goto cleanup;
    }


    double matrix_time =
        omp_get_wtime()
        - matrix_start;


    printf(
        "Harmonic matrix build time : %.6f s\n",
        matrix_time
    );


    /*
     * ========================================================
     * 11. LU
     * ========================================================
     */
    printf("\nFactorizing harmonic PEEC system ...\n");


    double LU_start =
        omp_get_wtime();


    if (harmonic_factorize(
            &harmonic_system,
            config.parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic LU factorization failed.\n"
        );

        goto cleanup;
    }


    double LU_time =
        omp_get_wtime()
        - LU_start;


    printf(
        "Harmonic LU time : %.6f s\n",
        LU_time
    );


    /*
     * ========================================================
     * 12. SOLUTION
     * ========================================================
     */
    HarmonicExcitation excitation;


    excitation.n_edges =
        mesh.n_edges;

    excitation.n_nodes =
        mesh.n_nodes;

    excitation.edge_voltage =
        U_inc;


    /*
     * При электромагнитном рассеянии:
     *
     *     I_s = 0.
     */
    excitation.node_current =
        NULL;


    printf("\nSolving harmonic PEEC system ...\n");


    double solve_start =
        omp_get_wtime();


    if (harmonic_solve(
            &harmonic_system,
            &P,
            &excitation,
            &solution,
            config.parallel_threads
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic solution failed.\n"
        );

        goto cleanup;
    }


    double solve_time =
        omp_get_wtime()
        - solve_start;


    /*
     * ========================================================
     * 13. RESIDUAL
     * ========================================================
     */
    HarmonicResidual residual;


    if (harmonic_compute_residual(
            &harmonic_system,
            &P,
            &excitation,
            &solution,
            config.parallel_threads,
            &residual
        ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: residual calculation failed.\n"
        );

        goto cleanup;
    }


    double maximum_current =
        complex_array_max_abs(
            solution.edge_current,
            solution.n_edges
        );


    double maximum_voltage =
        complex_array_max_abs(
            solution.node_voltage,
            solution.n_nodes
        );


    printf("\n");
    printf("============================================================\n");
    printf("HARMONIC PEEC SOLUTION\n");
    printf("============================================================\n");

    printf(
        "Maximum |I_e|      : %.9e A\n",
        maximum_current
    );

    printf(
        "Maximum |V_j|      : %.9e V\n",
        maximum_voltage
    );

    printf(
        "Backward error     : %.9e\n",
        residual.backward_error
    );

    printf(
        "Edge back. error   : %.9e\n",
        residual.edge_backward_error
    );

    printf(
        "Node back. error   : %.9e\n",
        residual.node_backward_error
    );

    printf(
        "Solve time         : %.6f s\n",
        solve_time
    );

    printf("============================================================\n");


    if (!isfinite(residual.backward_error) ||
        residual.backward_error > 1.0e-10)
    {
        fprintf(
            stderr,
            "ERROR: PEEC solution has a large backward error.\n"
        );

        goto cleanup;
    }


    /*
     * ========================================================
     * 14. RCS
     * ========================================================
     */
    RcsOptions rcs_options;


    rcs_options_default(
        &rcs_options,
        config.parallel_threads
    );


    rcs_options.integration =
        config.rcs_integration;

    rcs_options.surface_order =
        config.rcs_surface_order;


    /*
     * ========================================================
     * 15. OUTPUT DIRECTORY
     * ========================================================
     */
    if (ensure_results_directory() != 0)
    {
        goto cleanup;
    }


    /*
     * ========================================================
     * 16. RCS SWEEP
     * ========================================================
     */
    if (save_rcs_sweep(
            &config,
            &mesh,
            &dual,
            &field,
            &solution,
            &rcs_options
        ) != 0)
    {
        goto cleanup;
    }


    exit_code =
        EXIT_SUCCESS;


cleanup:

    free(
        U_inc
    );

    harmonic_solution_free(
        &solution
    );

    harmonic_system_free(
        &harmonic_system
    );

    resistance_free(
        &R
    );

    matrix_free(
        &PAT
    );

    potential_free(
        &P
    );

    inductance_free(
        &L
    );

    incidence_free(
        &incidence
    );

    dual_mesh_free(
        &dual
    );

    mesh_free(
        &mesh
    );


    if (exit_code == EXIT_SUCCESS)
    {
        printf("\n");
        printf("RCS calculation completed successfully.\n");
    }


    return exit_code;
}