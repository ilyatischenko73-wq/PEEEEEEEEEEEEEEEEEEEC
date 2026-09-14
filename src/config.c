#include "config.h"
#include "quadrature.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


/*
 * ============================================================
 * ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ПАРСИНГА
 * ============================================================
 */
struct PeecConfigString {
    struct PeecConfigString *next;
    char text[];
};

static int config_copy_string(PeecConfig *config, const char **field, const char *text)
{
    size_t length = strlen(text);
    struct PeecConfigString *copy = malloc(sizeof(*copy) + length + 1);
    if (copy == NULL) {
        fprintf(stderr, "ERROR: cannot allocate configuration string.\n");
        return -1;
    }
    memcpy(copy->text, text, length + 1);
    copy->next = config->owned_strings;
    config->owned_strings = copy;
    *field = copy->text;
    return 0;
}

void config_free(PeecConfig *config)
{
    if (config == NULL) return;
    struct PeecConfigString *p = config->owned_strings;
    while (p != NULL) {
        struct PeecConfigString *next = p->next;
        free(p);
        p = next;
    }
    memset(config, 0, sizeof(*config));
}

static int config_parse_double(
    const char *text,
    double *value
)
{
    if (text == NULL ||
        value == NULL)
    {
        return -1;
    }

    errno = 0;

    char *end = NULL;

    double parsed =
    strtod(
        text,
        &end
    );

    if (errno != 0 ||
        end == text ||
        *end != '\0' ||
        !isfinite(parsed))
    {
        return -1;
    }

    *value =
    parsed;

    return 0;
}


static int config_parse_int(
    const char *text,
    int *value
)
{
    if (text == NULL ||
        value == NULL)
    {
        return -1;
    }

    errno = 0;

    char *end = NULL;

    long parsed =
    strtol(
        text,
        &end,
        10
    );

    if (errno != 0 ||
        end == text ||
        *end != '\0')
    {
        return -1;
    }

    if (parsed < -2147483647L - 1L ||
        parsed > 2147483647L)
    {
        return -1;
    }

    *value =
    (int)parsed;

    return 0;
}


static int config_parse_xyz(
    int argc,
    char **argv,
    int *index,
    double xyz[3]
)
{
    if (argv == NULL ||
        index == NULL ||
        xyz == NULL)
    {
        return -1;
    }

    if (*index + 3 >= argc)
    {
        return -1;
    }

    for (int c = 0; c < 3; ++c)
    {
        if (config_parse_double(
            argv[*index + 1 + c],
            &xyz[c]
        ) != 0)
        {
            return -1;
        }
    }

    *index +=
    3;

    return 0;
}


/*
 * ============================================================
 * DEFAULTS
 * ============================================================
 */
void config_set_defaults(
    PeecConfig *config
)
{
    if (config == NULL)
    {
        return;
    }

    memset(
        config,
        0,
        sizeof(*config)
    );

    /*
     * --------------------------------------------------------
     * COMMON
     * --------------------------------------------------------
     */
    config->task =
    PEEC_TASK_MESH_INFO;

    config->physics =
    PEEC_PHYSICS_QUASISTATIC;

    config->mesh_file =
    NULL;

    config->parallel_threads =
    1;


    /*
     * --------------------------------------------------------
     * SCATTERING
     * --------------------------------------------------------
     */
    config->frequency =
    300.0e6;

    config->field_amplitude =
    1.0;

    config->theta_deg =
    90.0;

    config->phi_deg =
    45.0;

    config->phase_rad =
    0.0;

    config->polarization =
    PEEC_POLARIZATION_HORIZONTAL;


    /*
     * --------------------------------------------------------
     * RCS
     * --------------------------------------------------------
     */
    config->ka =
    1.0;

    config->characteristic_length =
    1.0;

    config->rcs_use_dual =
    1;

    config->rcs_order =
    6;

    config->rcs_observation_theta =
    90.0;

    config->rcs_phi_start =
    0.0;

    config->rcs_phi_end =
    360.0;

    config->rcs_phi_step =
    1.0;

    config->output_file =
    "results/rcs.csv";


    /*
     * --------------------------------------------------------
     * OUTPUT / CACHE
     * --------------------------------------------------------
     */
    config->save_matrices =
    0;

    config->matrix_directory =
    "cache";

    config->write_vtk =
    0;

    config->vtk_directory =
    "results/vtk";

    config->vtk_every =
    1;


    /*
     * --------------------------------------------------------
     * LIGHTNING
     * --------------------------------------------------------
     */
    config->lightning_case =
    PEEC_LIGHTNING_CLOSED_BODY;

    config->excitation =
    PEEC_EXCITATION_LIGHTNING_CURRENT;

    config->closed_mesh_file =
    NULL;

    config->open_mesh_file =
    NULL;

    config->aperture_map_file =
    NULL;

    config->slot_cells_file =
    NULL;

    config->slot_width =
    0.0;

    config->slot_wall_span =
    0.0;

    config->slot_wall_thickness =
    0.0;

    config->slot_epsilon_r =
    1.0;

    config->slot_mu_r =
    1.0;

    /*
     * Эти значения являются только безопасными стартовыми
     * значениями конфигурации.
     *
     * Конкретный расчет молнии должен явно задавать dt/t_end.
     */
    config->dt =
    0.0;

    config->t_end =
    0.0;

    config->time_order =
    2;

    /*
     * Те же параметры double-exponential waveform,
     * которые используются в lightning_source_default().
     */
    config->lightning_K =
    102900.0;

    config->lightning_alpha =
    1500.0;

    config->lightning_beta =
    1.0e6;

    config->lightning_delay =
    0.0;

    /*
     * Падающий double-exponential EM pulse.
     *
     * По умолчанию амплитуда 1 V/m.
     */
    config->pulse_K =
    1.0;

    config->pulse_alpha =
    1500.0;

    config->pulse_beta =
    1.0e6;

    config->strike_xyz_set =
    0;

    config->return_xyz_set =
    0;

    config->strike_node =
    -1;

    config->return_node =
    -1;

    config->strike_node_set =
    0;

    config->return_node_set =
    0;


    /*
     * --------------------------------------------------------
     * INTERNAL SHUNT
     * --------------------------------------------------------
     */
    config->use_internal_shunt =
    0;

    config->shunt_resistance =
    1000.0;

    config->shunt_a_xyz_set =
    0;

    config->shunt_b_xyz_set =
    0;

    config->shunt_a_node =
    -1;

    config->shunt_b_node =
    -1;

    config->shunt_a_node_set =
    0;

    config->shunt_b_node_set =
    0;
}


/*
 * ============================================================
 * HELP
 * ============================================================
 */
void config_print_help(
    const char *program_name
)
{
    printf("  --config FILE       read a case file (PeecSolver)\n");
    printf("  Aliases: --field-amplitude = --E; --characteristic-length = --a\n");
    printf("  --rcs-use-dual = --rcs-dual; --rcs-observation-theta = --obs-theta\n");
    printf("  --rcs-phi-start/end/step = --phi-start/end/step\n");

    if (program_name == NULL)
    {
        program_name =
        "PeecSolver";
    }

    printf("\n");
    printf("PEEC / Partial Equivalent Circuit solver\n");
    printf("\n");

    printf("Usage:\n");
    printf("\n");
    printf("  %s --task TASK [options]\n", program_name);
    printf("\n");


    printf("Tasks:\n");
    printf("\n");

    printf("  --task mesh-info\n");
    printf("      Read mesh and print topology information.\n");
    printf("\n");

    printf("  --task scattering\n");
    printf("      Harmonic scattering problem.\n");
    printf("\n");

    printf("  --task rcs\n");
    printf("      Harmonic scattering followed by RCS calculation.\n");
    printf("\n");

    printf("  --task lightning\n");
    printf("      Transient lightning problem.\n");
    printf("\n");


    printf("Common options:\n");
    printf("\n");

    printf("  --mesh FILE\n");
    printf("      Main mesh for mesh-info, scattering and RCS.\n");
    printf("\n");

    printf("  --physics quasistatic|retarded\n");
    printf("      PEEC interaction model.\n");
    printf("      Default: quasistatic.\n");
    printf("\n");

    printf("  --parallel N\n");
    printf("      Number of OpenMP threads.\n");
    printf("      Default: 1.\n");
    printf("\n");


    printf("Scattering / RCS options:\n");
    printf("\n");

    printf("  --frequency F\n");
    printf("      Frequency [Hz].\n");
    printf("      Default: 300e6.\n");
    printf("\n");

    printf("  --E VALUE\n");
    printf("      Incident electric-field amplitude [V/m].\n");
    printf("      Default: 1.\n");
    printf("\n");

    printf("  --theta DEG\n");
    printf("      Incident-wave theta angle [deg].\n");
    printf("\n");

    printf("  --phi DEG\n");
    printf("      Incident-wave phi angle [deg].\n");
    printf("\n");

    printf("  --phase RAD\n");
    printf("      Incident-wave phasor phase [rad].\n");
    printf("\n");

    printf("  --polarization vertical|horizontal\n");
    printf("      Incident-wave polarization.\n");
    printf("\n");


    printf("RCS options:\n");
    printf("\n");

    printf("  --ka VALUE\n");
    printf("      Electrical size k*a.\n");
    printf("\n");

    printf("  --a VALUE\n");
    printf("      Characteristic length a [m].\n");
    printf("\n");

    printf("  --rcs-dual\n");
    printf("      Use dual-support current integration.\n");
    printf("\n");

    printf("  --rcs-no-dual\n");
    printf("      Disable dual-support RCS integration.\n");
    printf("\n");

    printf("  --rcs-order N\n");
    printf("      Quadrature order for RCS integration.\n");
    printf("      Default: 6.\n");
    printf("\n");

    printf("  --output FILE\n");
    printf("      Output CSV file.\n");
    printf("\n");

    printf("  --obs-theta DEG\n");
    printf("      RCS observation theta [deg].\n");
    printf("\n");

    printf("  --phi-start DEG\n");
    printf("      RCS sweep start phi [deg].\n");
    printf("\n");

    printf("  --phi-end DEG\n");
    printf("      RCS sweep end phi [deg].\n");
    printf("\n");

    printf("  --phi-step DEG\n");
    printf("      RCS sweep step phi [deg].\n");
    printf("\n");


    printf("Output / cache options:\n");
    printf("\n");

    printf("  --save-matrices\n");
    printf("      Save L, A, P, PAT and R in Matrix Market format.\n");
    printf("\n");

    printf("  --matrix-dir DIR\n");
    printf("      Matrix output directory. Default: cache.\n");
    printf("\n");

    printf("  --vtk\n");
    printf("      Write ParaView Legacy VTK files.\n");
    printf("\n");

    printf("  --vtk-dir DIR\n");
    printf("      ParaView output directory. Default: results/vtk.\n");
    printf("\n");

    printf("  --vtk-every N\n");
    printf("      Transient VTK output stride. Default: 1.\n");
    printf("\n");


    printf("Lightning options:\n");
    printf("\n");

    printf("  --excitation lightning-current|incident-pulse\n");
    printf("      lightning-current : node current source between strike/return.\n");
    printf("      incident-pulse     : incident double-exponential plane-wave pulse.\n");
    printf("      Default: lightning-current.\n");
    printf("\n");

    printf("  --lightning-case closed|two-stage|slot-cells\n");
    printf("      closed     : one fully closed body.\n");
    printf("      two-stage  : closed Stage 1 + open Stage 2.\n");
    printf("      slot-cells : one open PEEC mesh + distributed LC aperture cells.\n");
    printf("\n");

    printf("  --closed-mesh FILE\n");
    printf("      Closed-body mesh used by Stage 1.\n");
    printf("\n");

    printf("  --open-mesh FILE\n");
    printf("      Open-body mesh used by Stage 2.\n");
    printf("\n");

    printf("  --aperture-map FILE\n");
    printf("      Stage-1 -> Stage-2 aperture mapping file.\n");
    printf("\n");

    printf("  --slot-cells FILE\n");
    printf("      Distributed slot-cell topology file.\n");
    printf("\n");

    printf("  --slot-width VALUE\n");
    printf("      Physical slot width w [m].\n");
    printf("\n");

    printf("  --slot-wall-span VALUE\n");
    printf("      Wall span b used by the coplanar-strip model [m].\n");
    printf("\n");

    printf("  --slot-wall-thickness VALUE\n");
    printf("      Conducting wall thickness t [m]. Default: 0.\n");
    printf("\n");

    printf("  --slot-epsilon-r VALUE\n");
    printf("      Relative permittivity of slot medium. Default: 1.\n");
    printf("\n");

    printf("  --slot-mu-r VALUE\n");
    printf("      Relative permeability of slot medium. Default: 1.\n");
    printf("\n");

    printf("  --dt VALUE\n");
    printf("      Time step [s].\n");
    printf("\n");

    printf("  --t-end VALUE\n");
    printf("      Final simulation time [s].\n");
    printf("\n");

    printf("  --time-order 1|2\n");
    printf("      1 : Backward Euler, first order.\n");
    printf("      2 : trapezoidal rule, second order.\n");
    printf("      Default: 2.\n");
    printf("\n");

    printf("  --lightning-K VALUE\n");
    printf("      Double-exponential amplitude K [A].\n");
    printf("\n");

    printf("  --lightning-alpha VALUE\n");
    printf("      Slow exponential coefficient alpha [1/s].\n");
    printf("\n");

    printf("  --lightning-beta VALUE\n");
    printf("      Fast exponential coefficient beta [1/s].\n");
    printf("\n");

    printf("  --lightning-delay VALUE\n");
    printf("      Lightning delay [s].\n");
    printf("\n");

    printf("  --pulse-K VALUE\n");
    printf("      Incident-pulse amplitude coefficient [V/m]. Default: 1.\n");
    printf("\n");

    printf("  --pulse-alpha VALUE\n");
    printf("      Incident-pulse slow exponential coefficient [1/s].\n");
    printf("      Default: 1500.\n");
    printf("\n");

    printf("  --pulse-beta VALUE\n");
    printf("      Incident-pulse fast exponential coefficient [1/s].\n");
    printf("      Default: 1e6.\n");
    printf("\n");

    printf("  For incident-pulse, --theta, --phi and --polarization\n");
    printf("  define propagation direction and field polarization.\n");
    printf("\n");

    printf("  --strike X Y Z\n");
    printf("      Strike-point coordinates [m].\n");
    printf("\n");

    printf("  --return X Y Z\n");
    printf("      Return-point coordinates [m].\n");
    printf("\n");

    printf("  --strike-node N\n");
    printf("      Strike node by internal solver index, 0-based.\n");
    printf("\n");

    printf("  --return-node N\n");
    printf("      Return node by internal solver index, 0-based.\n");
    printf("\n");

    printf("  --shunt-R VALUE\n");
    printf("      Enable Stage-2 resistive shunt and set R [Ohm].\n");
    printf("\n");

    printf("  --shunt-a X Y Z\n");
    printf("      Stage-2 shunt terminal A coordinates [m].\n");
    printf("\n");

    printf("  --shunt-b X Y Z\n");
    printf("      Stage-2 shunt terminal B coordinates [m].\n");
    printf("\n");

    printf("  --shunt-a-node N\n");
    printf("      Stage-2 shunt terminal A by internal 0-based node index.\n");
    printf("\n");

    printf("  --shunt-b-node N\n");
    printf("      Stage-2 shunt terminal B by internal 0-based node index.\n");
    printf("\n");


    printf("  -h, --help\n");
    printf("      Show this help.\n");
    printf("\n");
}


/*
 * ============================================================
 * FINAL VALIDATION
 * ============================================================
 */
static int config_validate(
    PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    if (config->parallel_threads < 1)
    {
        fprintf(
            stderr,
            "ERROR: --parallel must be >= 1.\n"
        );

        return -1;
    }


    /*
     * --------------------------------------------------------
     * MESH INFO / SCATTERING / RCS
     * --------------------------------------------------------
     */
    if (config->task == PEEC_TASK_MESH_INFO ||
        config->task == PEEC_TASK_SCATTERING ||
        config->task == PEEC_TASK_RCS)
    {
        if (config->mesh_file == NULL)
        {
            fprintf(
                stderr,
                "ERROR: this task requires --mesh FILE.\n"
            );

            return -1;
        }
    }


    /*
     * --------------------------------------------------------
     * SCATTERING
     * --------------------------------------------------------
     */
    if (config->task == PEEC_TASK_SCATTERING)
    {
        if (!isfinite(config->frequency) ||
            config->frequency <= 0.0)
        {
            fprintf(
                stderr,
                "ERROR: --frequency must be > 0.\n"
            );

            return -1;
        }

        if (!isfinite(config->field_amplitude) ||
            config->field_amplitude < 0.0)
        {
            fprintf(
                stderr,
                "ERROR: --E must be >= 0.\n"
            );

            return -1;
        }
    }


    /*
     * --------------------------------------------------------
     * RCS
     * --------------------------------------------------------
     */
    if (config->task == PEEC_TASK_RCS)
    {
        if (!isfinite(config->ka) ||
            config->ka <= 0.0)
        {
            fprintf(
                stderr,
                "ERROR: --ka must be > 0.\n"
            );

            return -1;
        }

        if (!isfinite(config->characteristic_length) ||
            config->characteristic_length <= 0.0)
        {
            fprintf(
                stderr,
                "ERROR: --a must be > 0.\n"
            );

            return -1;
        }

        if (!gauss_legendre_order_supported(config->rcs_order))
        {
            fprintf(
                stderr,
                "ERROR: unsupported --rcs-order; use 2,3,4,5,6,8,10,12,16.\n"
            );

            return -1;
        }

        if (config->output_file == NULL)
        {
            fprintf(
                stderr,
                "ERROR: RCS requires an output file.\n"
            );

            return -1;
        }

        if (!isfinite(config->rcs_phi_start) ||
            !isfinite(config->rcs_phi_end) ||
            config->rcs_phi_end < config->rcs_phi_start ||
            config->rcs_phi_end - config->rcs_phi_start > 1e6 * config->rcs_phi_step ||
            config->rcs_phi_start + config->rcs_phi_step == config->rcs_phi_start ||
            !isfinite(config->rcs_phi_step) ||
            config->rcs_phi_step <= 0.0)
        {
            fprintf(
                stderr,
                "ERROR: invalid RCS angle range/step (maximum 1000001 samples).\n"
            );

            return -1;
        }
    }

    if (config->vtk_every < 1)
    {
        fprintf(
            stderr,
            "ERROR: --vtk-every must be >= 1.\n"
        );

        return -1;
    }


    /*
     * --------------------------------------------------------
     * LIGHTNING
     * --------------------------------------------------------
     */
    if (config->task == PEEC_TASK_LIGHTNING)
    {
        /*
         * Mesh selection depends on the lightning model.
         *
         * closed / two-stage:
         *
         *     --mesh may be used as shorthand for --closed-mesh.
         *
         * slot-cells:
         *
         *     --mesh may be used as shorthand for --open-mesh.
         */
        if (config->lightning_case ==
            PEEC_LIGHTNING_SLOT_CELLS)
        {
            if (config->open_mesh_file == NULL &&
                config->mesh_file != NULL)
            {
                config->open_mesh_file =
                config->mesh_file;
            }

            if (config->open_mesh_file == NULL)
            {
                fprintf(
                    stderr,
                    "ERROR: slot-cells lightning requires --open-mesh FILE.\n"
                );

                return -1;
            }
        }
        else
        {
            if (config->closed_mesh_file == NULL &&
                config->mesh_file != NULL)
            {
                config->closed_mesh_file =
                config->mesh_file;
            }

            if (config->closed_mesh_file == NULL)
            {
                fprintf(
                    stderr,
                    "ERROR: lightning requires --closed-mesh FILE.\n"
                );

                return -1;
            }
        }

        if (!isfinite(config->dt) ||
            config->dt <= 0.0)
        {
            fprintf(
                stderr,
                "ERROR: lightning requires --dt > 0.\n"
            );

            return -1;
        }

        if (!isfinite(config->t_end) ||
            config->t_end <= 0.0)
        {
            fprintf(
                stderr,
                "ERROR: lightning requires --t-end > 0.\n"
            );

            return -1;
        }

        if (config->t_end < config->dt)
        {
            fprintf(
                stderr,
                "ERROR: --t-end must be >= --dt.\n"
            );

            return -1;
        }

        if (config->time_order != 1 &&
            config->time_order != 2)
        {
            fprintf(
                stderr,
                "ERROR: --time-order must be 1 or 2.\n"
            );

            return -1;
        }

        if (config->excitation ==
            PEEC_EXCITATION_LIGHTNING_CURRENT)
        {
            if (!isfinite(config->lightning_K) ||
                !isfinite(config->lightning_alpha) ||
                !isfinite(config->lightning_beta) ||
                !isfinite(config->lightning_delay) ||
                config->lightning_alpha <= 0.0 ||
                config->lightning_beta <= 0.0 ||
                config->lightning_delay < 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid lightning waveform parameters.\n"
                );

                return -1;
            }

            int lightning_by_xyz =
            config->strike_xyz_set
            && config->return_xyz_set;

            int lightning_by_node =
            config->strike_node_set
            && config->return_node_set;

            if (lightning_by_xyz &&
                lightning_by_node)
            {
                fprintf(
                    stderr,
                    "ERROR: choose either --strike/--return coordinates "
                    "or --strike-node/--return-node, not both.\n"
                );

                return -1;
            }

            if (!lightning_by_xyz &&
                !lightning_by_node)
            {
                fprintf(
                    stderr,
                    "ERROR: lightning-current requires either "
                    "--strike X Y Z + --return X Y Z, or "
                    "--strike-node N + --return-node N.\n"
                );

                return -1;
            }
        }
        else if (config->excitation ==
            PEEC_EXCITATION_INCIDENT_PULSE)
        {
            if (!isfinite(config->pulse_K) ||
                !isfinite(config->pulse_alpha) ||
                !isfinite(config->pulse_beta) ||
                config->pulse_alpha <= 0.0 ||
                config->pulse_beta <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid incident-pulse parameters.\n"
                );

                return -1;
            }
        }
        else
        {
            fprintf(
                stderr,
                "ERROR: invalid transient excitation.\n"
            );

            return -1;
        }

        if (config->lightning_case ==
            PEEC_LIGHTNING_TWO_STAGE)
        {
            if (config->open_mesh_file == NULL)
            {
                fprintf(
                    stderr,
                    "ERROR: two-stage lightning requires --open-mesh FILE.\n"
                );

                return -1;
            }

            if (config->aperture_map_file == NULL)
            {
                fprintf(
                    stderr,
                    "ERROR: two-stage lightning requires --aperture-map FILE.\n"
                );

                return -1;
            }
        }
        else if (config->lightning_case ==
            PEEC_LIGHTNING_SLOT_CELLS)
        {
            if (config->slot_cells_file == NULL)
            {
                fprintf(
                    stderr,
                    "ERROR: slot-cells lightning requires --slot-cells FILE.\n"
                );

                return -1;
            }

            if (!isfinite(config->slot_width) ||
                config->slot_width <= 0.0 ||
                !isfinite(config->slot_wall_span) ||
                config->slot_wall_span <= 0.0 ||
                !isfinite(config->slot_wall_thickness) ||
                config->slot_wall_thickness < 0.0 ||
                !isfinite(config->slot_epsilon_r) ||
                config->slot_epsilon_r <= 0.0 ||
                !isfinite(config->slot_mu_r) ||
                config->slot_mu_r <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid slot-cell geometry/material parameters.\n"
                );

                return -1;
            }
        }
        else
        {
            /*
             * Closed body has no internal shunt in the current driver.
             */
            config->use_internal_shunt =
            0;
        }

        /*
         * An internal shunt is supported in:
         *
         *     two-stage
         *     slot-cells.
         */
        if (config->use_internal_shunt)
        {
            if (config->lightning_case ==
                PEEC_LIGHTNING_CLOSED_BODY)
            {
                fprintf(
                    stderr,
                    "ERROR: internal shunt is not enabled for closed mode.\n"
                );

                return -1;
            }

            if (!isfinite(config->shunt_resistance) ||
                config->shunt_resistance < 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: --shunt-R must be >= 0 Ohm.\n"
                );

                return -1;
            }

            int shunt_by_xyz =
            config->shunt_a_xyz_set
            && config->shunt_b_xyz_set;

            int shunt_by_node =
            config->shunt_a_node_set
            && config->shunt_b_node_set;

            if (shunt_by_xyz &&
                shunt_by_node)
            {
                fprintf(
                    stderr,
                    "ERROR: choose either shunt coordinates or "
                    "--shunt-a-node/--shunt-b-node, not both.\n"
                );

                return -1;
            }

            if (!shunt_by_xyz &&
                !shunt_by_node)
            {
                fprintf(
                    stderr,
                    "ERROR: enabled shunt requires either "
                    "--shunt-a X Y Z + --shunt-b X Y Z, or "
                    "--shunt-a-node N + --shunt-b-node N.\n"
                );

                return -1;
            }
        }
    }

    return 0;
}


/*
 * ============================================================
 * CLI PARSER
 * ============================================================
 */
int config_parse_cli(
    PeecConfig *config,
    int argc,
    char **argv
)
{
    if (config == NULL ||
        argv == NULL)
    {
        return -1;
    }

    for (int i = 1; i < argc; ++i)
    {
        const char *arg =
        argv[i];


        /*
         * ----------------------------------------------------
         * HELP
         * ----------------------------------------------------
         */
        if (strcmp(arg, "-h") == 0 ||
            strcmp(arg, "--help") == 0)
        {
            config_print_help(
                argv[0]
            );

            return 1;
        }


        /*
         * ----------------------------------------------------
         * TASK
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--task") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --task requires a value.\n"
                );

                return -1;
            }

            const char *value =
            argv[++i];

            if (strcmp(value, "mesh-info") == 0)
            {
                config->task =
                PEEC_TASK_MESH_INFO;
            }
            else if (strcmp(value, "scattering") == 0)
            {
                config->task =
                PEEC_TASK_SCATTERING;
            }
            else if (strcmp(value, "rcs") == 0)
            {
                config->task =
                PEEC_TASK_RCS;
            }
            else if (strcmp(value, "lightning") == 0)
            {
                config->task =
                PEEC_TASK_LIGHTNING;
            }
            else
            {
                fprintf(
                    stderr,
                    "ERROR: unknown task: %s\n",
                    value
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * PHYSICS
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--physics") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --physics requires a value.\n"
                );

                return -1;
            }

            const char *value =
            argv[++i];

            if (strcmp(value, "quasistatic") == 0)
            {
                config->physics =
                PEEC_PHYSICS_QUASISTATIC;
            }
            else if (strcmp(value, "retarded") == 0)
            {
                config->physics =
                PEEC_PHYSICS_RETARDED;
            }
            else
            {
                fprintf(
                    stderr,
                    "ERROR: unknown physics mode: %s\n",
                    value
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * COMMON MESH
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--mesh") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --mesh requires a filename.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->mesh_file, argv[++i]) != 0) return -1;
        }


        /*
         * ----------------------------------------------------
         * PARALLEL
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--parallel") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->parallel_threads
                ) != 0 ||
                config->parallel_threads < 1)
            {
                fprintf(
                    stderr,
                    "ERROR: --parallel must be an integer >= 1.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * SCATTERING / RCS FIELD
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--frequency") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->frequency
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --frequency value.\n"
                );

                return -1;
            }
        }

        else if ((strcmp(arg, "--E") == 0 || strcmp(arg, "--field-amplitude") == 0))
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->field_amplitude
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --E value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--theta") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->theta_deg
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --theta value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--phi") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->phi_deg
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phi value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--phase") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->phase_rad
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phase value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--polarization") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --polarization requires a value.\n"
                );

                return -1;
            }

            const char *value =
            argv[++i];

            if (strcmp(value, "vertical") == 0)
            {
                config->polarization =
                PEEC_POLARIZATION_VERTICAL;
            }
            else if (strcmp(value, "horizontal") == 0)
            {
                config->polarization =
                PEEC_POLARIZATION_HORIZONTAL;
            }
            else
            {
                fprintf(
                    stderr,
                    "ERROR: unknown polarization: %s\n",
                    value
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * RCS
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--ka") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->ka
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --ka value.\n"
                );

                return -1;
            }
        }

        else if ((strcmp(arg, "--a") == 0 || strcmp(arg, "--characteristic-length") == 0))
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->characteristic_length
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --a value.\n"
                );

                return -1;
            }
        }

        else if ((strcmp(arg, "--rcs-dual") == 0 || strcmp(arg, "--rcs-use-dual") == 0))
        {
            config->rcs_use_dual =
            1;
        }

        else if (strcmp(arg, "--rcs-no-dual") == 0)
        {
            config->rcs_use_dual =
            0;
        }

        else if (strcmp(arg, "--rcs-order") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->rcs_order
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --rcs-order value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--output") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --output requires a filename.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->output_file, argv[++i]) != 0) return -1;
        }


        /*
         * ----------------------------------------------------
         * RCS SWEEP
         * ----------------------------------------------------
         */
        else if ((strcmp(arg, "--obs-theta") == 0 || strcmp(arg, "--rcs-observation-theta") == 0))
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->rcs_observation_theta
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --obs-theta value.\n"
                );

                return -1;
            }
        }

        else if ((strcmp(arg, "--phi-start") == 0 || strcmp(arg, "--rcs-phi-start") == 0))
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->rcs_phi_start
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phi-start value.\n"
                );

                return -1;
            }
        }

        else if ((strcmp(arg, "--phi-end") == 0 || strcmp(arg, "--rcs-phi-end") == 0))
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->rcs_phi_end
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --phi-end value.\n"
                );

                return -1;
            }
        }

        else if ((strcmp(arg, "--phi-step") == 0 || strcmp(arg, "--rcs-phi-step") == 0))
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->rcs_phi_step
                ) != 0 ||
                config->rcs_phi_step <= 0.0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid RCS angle range/step (maximum 1000001 samples).\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * OUTPUT / CACHE
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--save-matrices") == 0)
        {
            config->save_matrices =
            1;
        }

        else if (strcmp(arg, "--matrix-dir") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --matrix-dir requires a directory.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->matrix_directory, argv[++i]) != 0) return -1;
        }

        else if (strcmp(arg, "--vtk") == 0)
        {
            config->write_vtk =
            1;
        }

        else if (strcmp(arg, "--vtk-dir") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --vtk-dir requires a directory.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->vtk_directory, argv[++i]) != 0) return -1;
        }

        else if (strcmp(arg, "--vtk-every") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->vtk_every
                ) != 0 ||
                config->vtk_every < 1)
            {
                fprintf(
                    stderr,
                    "ERROR: --vtk-every must be an integer >= 1.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * LIGHTNING CASE
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--excitation") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --excitation requires a value.\n"
                );

                return -1;
            }

            const char *value =
            argv[++i];

            if (strcmp(value, "lightning-current") == 0)
            {
                config->excitation =
                PEEC_EXCITATION_LIGHTNING_CURRENT;
            }
            else if (strcmp(value, "incident-pulse") == 0)
            {
                config->excitation =
                PEEC_EXCITATION_INCIDENT_PULSE;
            }
            else
            {
                fprintf(
                    stderr,
                    "ERROR: unknown excitation: %s\n",
                    value
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--lightning-case") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --lightning-case requires a value.\n"
                );

                return -1;
            }

            const char *value =
            argv[++i];

            if (strcmp(value, "closed") == 0)
            {
                config->lightning_case =
                PEEC_LIGHTNING_CLOSED_BODY;
            }
            else if (strcmp(value, "two-stage") == 0)
            {
                config->lightning_case =
                PEEC_LIGHTNING_TWO_STAGE;
            }
            else if (strcmp(value, "slot-cells") == 0)
            {
                config->lightning_case =
                PEEC_LIGHTNING_SLOT_CELLS;
            }
            else
            {
                fprintf(
                    stderr,
                    "ERROR: unknown lightning case: %s\n",
                    value
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--closed-mesh") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --closed-mesh requires a filename.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->closed_mesh_file, argv[++i]) != 0) return -1;
        }

        else if (strcmp(arg, "--open-mesh") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --open-mesh requires a filename.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->open_mesh_file, argv[++i]) != 0) return -1;
        }

        else if (strcmp(arg, "--aperture-map") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --aperture-map requires a filename.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->aperture_map_file, argv[++i]) != 0) return -1;
        }


        else if (strcmp(arg, "--slot-cells") == 0)
        {
            if (i + 1 >= argc)
            {
                fprintf(
                    stderr,
                    "ERROR: --slot-cells requires a filename.\n"
                );

                return -1;
            }

            if (config_copy_string(config, &config->slot_cells_file, argv[++i]) != 0) return -1;
        }

        else if (strcmp(arg, "--slot-width") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->slot_width
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --slot-width value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--slot-wall-span") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->slot_wall_span
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --slot-wall-span value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--slot-wall-thickness") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->slot_wall_thickness
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --slot-wall-thickness value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--slot-epsilon-r") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->slot_epsilon_r
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --slot-epsilon-r value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--slot-mu-r") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->slot_mu_r
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --slot-mu-r value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * LIGHTNING TIME GRID
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--dt") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->dt
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --dt value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--t-end") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->t_end
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --t-end value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--time-order") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->time_order
                ) != 0 ||
                (
                    config->time_order != 1 &&
                    config->time_order != 2
                ))
            {
                fprintf(
                    stderr,
                    "ERROR: --time-order must be 1 or 2.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * LIGHTNING WAVEFORM
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--lightning-K") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->lightning_K
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --lightning-K value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--lightning-alpha") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->lightning_alpha
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --lightning-alpha value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--lightning-beta") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->lightning_beta
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --lightning-beta value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--lightning-delay") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->lightning_delay
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --lightning-delay value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * INCIDENT PULSE
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--pulse-K") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->pulse_K
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --pulse-K value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--pulse-alpha") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->pulse_alpha
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --pulse-alpha value.\n"
                );

                return -1;
            }
        }

        else if (strcmp(arg, "--pulse-beta") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->pulse_beta
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --pulse-beta value.\n"
                );

                return -1;
            }
        }


        /*
         * ----------------------------------------------------
         * LIGHTNING TERMINALS
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--strike") == 0)
        {
            if (config_parse_xyz(
                argc,
                argv,
                &i,
                config->strike_xyz
            ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --strike requires X Y Z.\n"
                );

                return -1;
            }

            config->strike_xyz_set =
            1;
        }

        else if (strcmp(arg, "--return") == 0)
        {
            if (config_parse_xyz(
                argc,
                argv,
                &i,
                config->return_xyz
            ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --return requires X Y Z.\n"
                );

                return -1;
            }

            config->return_xyz_set =
            1;
        }

        else if (strcmp(arg, "--strike-node") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->strike_node
                ) != 0 ||
                config->strike_node < 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --strike-node must be an integer >= 0.\n"
                );

                return -1;
            }

            config->strike_node_set =
            1;
        }

        else if (strcmp(arg, "--return-node") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->return_node
                ) != 0 ||
                config->return_node < 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --return-node must be an integer >= 0.\n"
                );

                return -1;
            }

            config->return_node_set =
            1;
        }


        /*
         * ----------------------------------------------------
         * STAGE-2 SHUNT
         * ----------------------------------------------------
         */
        else if (strcmp(arg, "--shunt-R") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_double(
                    argv[++i],
                    &config->shunt_resistance
                ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: invalid --shunt-R value.\n"
                );

                return -1;
            }

            config->use_internal_shunt =
            1;
        }

        else if (strcmp(arg, "--shunt-a") == 0)
        {
            if (config_parse_xyz(
                argc,
                argv,
                &i,
                config->shunt_a_xyz
            ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --shunt-a requires X Y Z.\n"
                );

                return -1;
            }

            config->shunt_a_xyz_set =
            1;
        }

        else if (strcmp(arg, "--shunt-b") == 0)
        {
            if (config_parse_xyz(
                argc,
                argv,
                &i,
                config->shunt_b_xyz
            ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --shunt-b requires X Y Z.\n"
                );

                return -1;
            }

            config->shunt_b_xyz_set =
            1;
        }

        else if (strcmp(arg, "--shunt-a-node") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->shunt_a_node
                ) != 0 ||
                config->shunt_a_node < 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --shunt-a-node must be an integer >= 0.\n"
                );

                return -1;
            }

            config->shunt_a_node_set =
            1;
        }

        else if (strcmp(arg, "--shunt-b-node") == 0)
        {
            if (i + 1 >= argc ||
                config_parse_int(
                    argv[++i],
                    &config->shunt_b_node
                ) != 0 ||
                config->shunt_b_node < 0)
            {
                fprintf(
                    stderr,
                    "ERROR: --shunt-b-node must be an integer >= 0.\n"
                );

                return -1;
            }

            config->shunt_b_node_set =
            1;
        }


        /*
         * ----------------------------------------------------
         * UNKNOWN
         * ----------------------------------------------------
         */
        else
        {
            fprintf(
                stderr,
                "ERROR: unknown argument: %s\n",
                arg
            );

            return -1;
        }
    }

    return config_validate(
        config
    );
}
