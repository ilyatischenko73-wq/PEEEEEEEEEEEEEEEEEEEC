#include "config.h"

#include "mesh.h"
#include "topology.h"
#include "dual_mesh.h"

#include "incidence.h"
#include "inductance.h"
#include "potential.h"
#include "resistance.h"
#include "matrix.h"

#include "lu.h"
#include "harmonic.h"
#include "transient.h"

#include "incident_field.h"
#include "scattering.h"
#include "rcs.h"

#include "lightning.h"
#include "internal_circuit.h"
#include "aperture_coupling.h"
#include "slot_circuit.h"

#include "matrix_export.h"
#include "paraview.h"

#include <errno.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#endif


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


static const double C0 =
299792458.0;

static const double EPSILON0 =
8.8541878128e-12;

static const double MU0 =
1.2566370614359172953850573533118e-6;


/*
 * ============================================================
 * PEEC MODEL
 * ============================================================
 *
 * Один полностью собранный набор:
 *
 *     mesh
 *     dual
 *     A
 *     L
 *     P
 *     P A^T
 *     R.
 *
 * PAT сохраняется для harmonic solver.
 *
 * Transient solver использует новую ЧЭС-запись через:
 *
 *     F_j = 1/P_jj
 *     AS  = A S.
 *
 *
 * Для two-stage lightning создаются два независимых объекта:
 *
 *     closed_model
 *     open_model.
 */
typedef struct
{
    Mesh mesh;
    DualMesh dual;

    IncidenceMatrix incidence;

    InductanceMatrix L;
    PotentialMatrix P;
    Matrix PAT;
    ResistanceMatrix R;

    int assembled;

} PeecModel;


/*
 * ============================================================
 * CHARGE RECOVERY
 * ============================================================
 *
 * Решатель программы использует неизвестные:
 *
 *     I, V.
 *
 *
 * Для визуализации заряда восстанавливаем:
 *
 *     P Q = V
 *
 *     Q = P^{-1} V.
 *
 *
 * P^{-1} явно не вычисляется.
 * LU(P) строится один раз.
 */
typedef struct
{
    LUFactorization lu;
    size_t n;
    int factorized;

} ChargeRecovery;


/*
 * ============================================================
 * NAMES
 * ============================================================
 */
static const char *task_name(
    PeecTask task
)
{
    switch (task)
    {
        case PEEC_TASK_MESH_INFO:
            return "mesh-info";

        case PEEC_TASK_SCATTERING:
            return "scattering";

        case PEEC_TASK_RCS:
            return "rcs";

        case PEEC_TASK_LIGHTNING:
            return "lightning";

        case PEEC_TASK_NONE:
        default:
            return "none";
    }
}


static const char *physics_name(
    PeecPhysicsMode physics
)
{
    switch (physics)
    {
        case PEEC_PHYSICS_QUASISTATIC:
            return "quasistatic";

        case PEEC_PHYSICS_RETARDED:
            return "retarded";

        default:
            return "unknown";
    }
}


static TransientScheme transient_scheme_from_config(
    const PeecConfig *config
)
{
    if (config != NULL &&
        config->time_order == 1)
    {
        return TRANSIENT_BACKWARD_EULER;
    }

    return TRANSIENT_TRAPEZOIDAL;
}


static const char *excitation_name(
    PeecExcitation excitation
)
{
    switch (excitation)
    {
        case PEEC_EXCITATION_LIGHTNING_CURRENT:
            return "lightning-current";

        case PEEC_EXCITATION_INCIDENT_PULSE:
            return "incident-pulse";

        default:
            return "unknown";
    }
}


/*
 * ============================================================
 * CONSTANT-DT TIME GRID
 * ============================================================
 *
 * t_end должно быть целым числом шагов dt:
 *
 *     t_end = N * dt.
 *
 *
 * Нельзя использовать:
 *
 *     ceil(t_end / dt),
 *
 * потому что из-за двоичного представления double, например:
 *
 *     20e-6 / 1e-6
 *
 * может получиться:
 *
 *     20.000000000000004,
 *
 * и ceil() ошибочно даст 21.
 */
static int compute_constant_dt_steps(
    double t_end,
    double dt,
    size_t *n_steps
)
{
    if (n_steps == NULL ||
        !isfinite(t_end) ||
        !isfinite(dt) ||
        t_end <= 0.0 ||
        dt <= 0.0)
    {
        return -1;
    }

    double ratio =
    t_end / dt;

    long long rounded =
    llround(ratio);

    if (rounded < 1)
    {
        return -1;
    }

    double error =
    fabs(
        ratio
        - (double)rounded
    );

    double tolerance =
    1.0e-10
    * fmax(
        1.0,
        fabs(ratio)
    );

    if (error > tolerance)
    {
        fprintf(
            stderr,
            "ERROR: t_end must be an integer multiple of dt. "
            "t_end/dt = %.17g, nearest integer = %lld.\n",
            ratio,
            rounded
        );

        return -1;
    }

    *n_steps =
    (size_t)rounded;

    return 0;
}


/*
 * ============================================================
 * FILESYSTEM HELPERS
 * ============================================================
 */
static int peec_mkdir(
    const char *path
)
{
    #ifdef _WIN32
    return _mkdir(
        path
    );
    #else
    return mkdir(
        path,
        0755
    );
    #endif
}


static int is_path_separator(
    char c
)
{
    return
    c == '/'
    ||
    c == '\\';
}


static int ensure_directory(
    const char *directory
)
{
    if (directory == NULL ||
        directory[0] == '\0')
    {
        return -1;
    }

    char buffer[4096];

    size_t length =
    strlen(directory);

    if (length >= sizeof(buffer))
    {
        return -1;
    }

    memcpy(
        buffer,
        directory,
        length + 1
    );

    for (size_t i = 1; i <= length; ++i)
    {
        if (is_path_separator(buffer[i]) ||
            buffer[i] == '\0')
        {
            char saved =
            buffer[i];

            buffer[i] =
            '\0';

            if (buffer[0] != '\0' &&
                peec_mkdir(buffer) != 0 &&
                errno != EEXIST)
            {
                fprintf(
                    stderr,
                    "ERROR: cannot create directory '%s': %s\n",
                    buffer,
                    strerror(errno)
                );

                return -1;
            }

            buffer[i] =
            saved;
        }
    }

    return 0;
}


static int ensure_parent_directory(
    const char *filename
)
{
    if (filename == NULL)
    {
        return -1;
    }

    char buffer[4096];

    size_t length =
    strlen(filename);

    if (length >= sizeof(buffer))
    {
        return -1;
    }

    memcpy(
        buffer,
        filename,
        length + 1
    );

    char *slash_forward =
    strrchr(
        buffer,
        '/'
    );

    char *slash_backward =
    strrchr(
        buffer,
        '\\'
    );

    char *slash =
    slash_forward;

    if (slash_backward != NULL &&
        (
            slash == NULL ||
            slash_backward > slash
        ))
    {
        slash =
        slash_backward;
    }

    if (slash == NULL)
    {
        return 0;
    }

    *slash =
    '\0';

    if (buffer[0] == '\0')
    {
        return 0;
    }

    return ensure_directory(
        buffer
    );
}


static int make_subdirectory(
    char *buffer,
    size_t size,
    const char *base,
    const char *name
)
{
    if (buffer == NULL ||
        base == NULL ||
        name == NULL)
    {
        return -1;
    }

    int written =
    snprintf(
        buffer,
        size,
        "%s/%s",
        base,
        name
    );

    if (written < 0 ||
        (size_t)written >= size)
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * PEEC MODEL INIT / FREE
 * ============================================================
 */
static void peec_model_init(
    PeecModel *model
)
{
    if (model == NULL)
    {
        return;
    }

    memset(
        model,
        0,
        sizeof(*model)
    );

    mesh_init(
        &model->mesh
    );

    dual_mesh_init(
        &model->dual
    );

    incidence_init(
        &model->incidence
    );

    inductance_init(
        &model->L
    );

    potential_init(
        &model->P
    );

    matrix_init(
        &model->PAT
    );

    resistance_init(
        &model->R
    );
}


static void peec_model_free(
    PeecModel *model
)
{
    if (model == NULL)
    {
        return;
    }

    resistance_free(
        &model->R
    );

    matrix_free(
        &model->PAT
    );

    potential_free(
        &model->P
    );

    inductance_free(
        &model->L
    );

    incidence_free(
        &model->incidence
    );

    dual_mesh_free(
        &model->dual
    );

    mesh_free(
        &model->mesh
    );

    memset(
        model,
        0,
        sizeof(*model)
    );
}


/*
 * ============================================================
 * MESH ONLY
 * ============================================================
 */
static int peec_load_mesh_geometry(
    const char *mesh_file,
    PeecModel *model
)
{
    if (mesh_file == NULL ||
        model == NULL)
    {
        return -1;
    }

    printf("\nReading Gmsh mesh ...\n");

    if (mesh_load_gmsh(
        mesh_file,
        &model->mesh
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: mesh loading failed: %s\n",
            mesh_file
        );

        return -1;
    }

    printf("Building topology ...\n");

    if (mesh_build_topology(
        &model->mesh
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: topology construction failed.\n"
        );

        return -1;
    }

    mesh_print_info(
        &model->mesh
    );

    printf("\nBuilding dual mesh ...\n");

    if (dual_mesh_build(
        &model->mesh,
        &model->dual
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: dual mesh construction failed.\n"
        );

        return -1;
    }

    dual_mesh_print_info(
        &model->dual
    );

    return 0;
}


/*
 * ============================================================
 * FULL PEEC ASSEMBLY
 * ============================================================
 */
static int peec_assemble_model(
    const char *mesh_file,
    int parallel_threads,
    PeecModel *model
)
{
    if (mesh_file == NULL ||
        model == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (peec_load_mesh_geometry(
        mesh_file,
        model
    ) != 0)
    {
        return -1;
    }


    /*
     * --------------------------------------------------------
     * A
     * --------------------------------------------------------
     */
    printf("\nBuilding incidence matrix A ...\n");

    if (incidence_build(
        &model->mesh,
        &model->incidence
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: incidence matrix construction failed.\n"
        );

        return -1;
    }


    /*
     * --------------------------------------------------------
     * L
     * --------------------------------------------------------
     */
    InductanceOptions inductance_options;

    inductance_options_default(
        &inductance_options,
        parallel_threads
    );

    printf("\nComputing inductance matrix L ...\n");

    double L_start =
    omp_get_wtime();

    if (inductance_compute_matrix(
        &model->mesh,
        &model->dual,
        &inductance_options,
        &model->L
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: inductance matrix calculation failed.\n"
        );

        return -1;
    }

    double L_time =
    omp_get_wtime()
    - L_start;


    /*
     * --------------------------------------------------------
     * P
     * --------------------------------------------------------
     */
    PotentialOptions potential_options;

    potential_options_default(
        &potential_options,
        parallel_threads
    );

    printf("\nComputing potential matrix P ...\n");

    double P_start =
    omp_get_wtime();

    if (potential_compute_matrix(
        &model->mesh,
        &model->dual,
        &potential_options,
        &model->P
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: potential matrix calculation failed.\n"
        );

        return -1;
    }

    double P_time =
    omp_get_wtime()
    - P_start;


    /*
     * --------------------------------------------------------
     * P A^T
     * --------------------------------------------------------
     */
    printf("\nBuilding P A^T ...\n");

    double PAT_start =
    omp_get_wtime();

    if (incidence_build_PAT(
        &model->mesh,
        &model->P,
        &model->PAT,
        parallel_threads
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: P A^T construction failed.\n"
        );

        return -1;
    }

    double PAT_time =
    omp_get_wtime()
    - PAT_start;


    /*
     * --------------------------------------------------------
     * R
     * --------------------------------------------------------
     */
    ResistanceOptions resistance_options;

    resistance_options_default(
        &resistance_options,
        parallel_threads
    );

    /*
     * Основная текущая модель:
     *
     *     PEC.
     */
    resistance_options.model =
    PEEC_CONDUCTOR_PEC;

    printf("\nBuilding resistance matrix R ...\n");

    double R_start =
    omp_get_wtime();

    if (resistance_compute_matrix(
        &model->mesh,
        &model->dual,
        &resistance_options,
        &model->R
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: resistance matrix calculation failed.\n"
        );

        return -1;
    }

    double R_time =
    omp_get_wtime()
    - R_start;


    printf("\n");
    printf("PEEC matrix assembly summary\n");
    printf("------------------------------------------------------------\n");

    printf(
        "L time             : %.6f s\n",
        L_time
    );

    printf(
        "P time             : %.6f s\n",
        P_time
    );

    printf(
        "P A^T time         : %.6f s\n",
        PAT_time
    );

    printf(
        "R time             : %.6f s\n",
        R_time
    );

    printf(
        "Total              : %.6f s\n",
        L_time
        + P_time
        + PAT_time
        + R_time
    );

    printf(
        "Full system size   : %zu x %zu\n",
        model->mesh.n_edges + model->mesh.n_nodes,
        model->mesh.n_edges + model->mesh.n_nodes
    );

    model->assembled =
    1;

    return 0;
}


/*
 * ============================================================
 * MATRIX EXPORT
 * ============================================================
 */
static int peec_export_matrices_if_requested(
    const PeecConfig *config,
    const PeecModel *model,
    const char *directory
)
{
    if (config == NULL ||
        model == NULL ||
        directory == NULL)
    {
        return -1;
    }

    if (!config->save_matrices)
    {
        return 0;
    }

    return matrix_export_peec_set(
        directory,
        &model->mesh,
        &model->L,
        &model->P,
        &model->PAT,
        &model->R
    );
}


/*
 * ============================================================
 * CHARGE RECOVERY
 * ============================================================
 */
static void charge_recovery_init(
    ChargeRecovery *recovery
)
{
    if (recovery == NULL)
    {
        return;
    }

    memset(
        recovery,
        0,
        sizeof(*recovery)
    );

    lu_init(
        &recovery->lu
    );
}


static void charge_recovery_free(
    ChargeRecovery *recovery
)
{
    if (recovery == NULL)
    {
        return;
    }

    lu_free(
        &recovery->lu
    );

    memset(
        recovery,
        0,
        sizeof(*recovery)
    );
}


static int charge_recovery_factorize(
    ChargeRecovery *recovery,
    const PotentialMatrix *P,
    int parallel_threads
)
{
    if (recovery == NULL ||
        P == NULL ||
        P->data == NULL ||
        P->n == 0)
    {
        return -1;
    }

    lu_free(
        &recovery->lu
    );

    lu_init(
        &recovery->lu
    );

    LUOptions options;

    lu_options_default(
        &options,
        parallel_threads
    );

    if (lu_factorize(
        P->data,
        P->n,
        &options,
        &recovery->lu
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: LU(P) failed during charge recovery setup.\n"
        );

        return -1;
    }

    recovery->n =
    P->n;

    recovery->factorized =
    1;

    return 0;
}


static int charge_recovery_real(
    const ChargeRecovery *recovery,
    const double *V,
    double *Q
)
{
    if (recovery == NULL ||
        !recovery->factorized ||
        V == NULL ||
        Q == NULL)
    {
        return -1;
    }

    return lu_solve(
        &recovery->lu,
        V,
        Q
    );
}


static int charge_recovery_complex(
    const ChargeRecovery *recovery,
    const Complex *V,
    Complex *Q
)
{
    if (recovery == NULL ||
        !recovery->factorized ||
        V == NULL ||
        Q == NULL)
    {
        return -1;
    }

    size_t n =
    recovery->n;

    double *rhs =
    malloc(
        n * sizeof(double)
    );

    double *solution =
    malloc(
        n * sizeof(double)
    );

    if (rhs == NULL ||
        solution == NULL)
    {
        free(solution);
        free(rhs);

        return -1;
    }

    for (size_t j = 0; j < n; ++j)
    {
        rhs[j] =
        V[j].re;
    }

    if (charge_recovery_real(
        recovery,
        rhs,
        solution
    ) != 0)
    {
        free(solution);
        free(rhs);

        return -1;
    }

    for (size_t j = 0; j < n; ++j)
    {
        Q[j].re =
        solution[j];
    }

    for (size_t j = 0; j < n; ++j)
    {
        rhs[j] =
        V[j].im;
    }

    if (charge_recovery_real(
        recovery,
        rhs,
        solution
    ) != 0)
    {
        free(solution);
        free(rhs);

        return -1;
    }

    for (size_t j = 0; j < n; ++j)
    {
        Q[j].im =
        solution[j];
    }

    free(solution);
    free(rhs);

    return 0;
}


/*
 * ============================================================
 * COMPLEX VECTOR MAX
 * ============================================================
 */
static double complex_vector_max_abs(
    const Complex *x,
    size_t n
)
{
    if (x == NULL ||
        n == 0)
    {
        return 0.0;
    }

    double maximum =
    0.0;

    for (size_t i = 0; i < n; ++i)
    {
        double value =
        hypot(
            x[i].re,
            x[i].im
        );

        if (value > maximum)
        {
            maximum =
            value;
        }
    }

    return maximum;
}


/*
 * ============================================================
 * INCIDENT POLARIZATION
 * ============================================================
 */
static IncidentPolarization incident_polarization_from_config(
    const PeecConfig *config
)
{
    if (config != NULL &&
        config->polarization ==
        PEEC_POLARIZATION_VERTICAL)
    {
        return INCIDENT_POLARIZATION_VERTICAL;
    }

    return INCIDENT_POLARIZATION_HORIZONTAL;
}


static int configure_incident_pulse(
    const PeecConfig *config,
    IncidentField *field
)
{
    if (config == NULL ||
        field == NULL)
    {
        return -1;
    }

    if (incident_field_configure_pulse(
        field,
        incident_polarization_from_config(
            config
        ),
        config->pulse_K,
        config->pulse_alpha,
        config->pulse_beta,
        config->theta_deg,
        config->phi_deg
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot configure incident pulse.\n"
        );

        return -1;
    }

    incident_field_print_info(
        field
    );

    return 0;
}


/*
 * ============================================================
 * HARMONIC SOLVE
 * ============================================================
 */
static int solve_harmonic_scattering(
    const PeecConfig *config,
    const PeecModel *model,
    double frequency,
    IncidentField *field,
    HarmonicSolution *solution
)
{
    if (config == NULL ||
        model == NULL ||
        field == NULL ||
        solution == NULL ||
        frequency <= 0.0)
    {
        return -1;
    }

    if (incident_field_configure_harmonic(
        field,
        INCIDENT_WAVE_COMPLEX,
        incident_polarization_from_config(
            config
        ),
        config->field_amplitude,
        frequency,
        config->theta_deg,
        config->phi_deg,
        config->phase_rad
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot configure incident field.\n"
        );

        return -1;
    }

    incident_field_print_info(
        field
    );


    /*
     * --------------------------------------------------------
     * EDGE EXCITATION
     * --------------------------------------------------------
     */
    Complex *U_inc =
    malloc(
        model->mesh.n_edges
        * sizeof(Complex)
    );

    if (U_inc == NULL)
    {
        return -1;
    }

    ScatteringOptions scattering_options;

    scattering_options_default(
        &scattering_options,
        config->parallel_threads
    );

    scattering_options.integration =
    SCATTERING_EDGE_ANALYTIC_PLANE_WAVE;

    if (scattering_build_harmonic_excitation(
        &model->mesh,
        field,
        &scattering_options,
        U_inc
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic scattering excitation failed.\n"
        );

        free(
            U_inc
        );

        return -1;
    }


    /*
     * --------------------------------------------------------
     * SYSTEM
     * --------------------------------------------------------
     */
    HarmonicSystem harmonic_system;

    harmonic_system_init(
        &harmonic_system
    );

    if (harmonic_build_system(
        &model->L,
        &model->R,
        &model->incidence,
        &model->PAT,
        frequency,
        config->parallel_threads,
        &harmonic_system
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic PEEC system construction failed.\n"
        );

        harmonic_system_free(
            &harmonic_system
        );

        free(
            U_inc
        );

        return -1;
    }

    if (harmonic_factorize(
        &harmonic_system,
        config->parallel_threads
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic factorization failed.\n"
        );

        harmonic_system_free(
            &harmonic_system
        );

        free(
            U_inc
        );

        return -1;
    }


    /*
     * --------------------------------------------------------
     * EXCITATION
     * --------------------------------------------------------
     */
    HarmonicExcitation excitation;

    excitation.n_edges =
    model->mesh.n_edges;

    excitation.n_nodes =
    model->mesh.n_nodes;

    excitation.edge_voltage =
    U_inc;

    excitation.node_current =
    NULL;


    /*
     * --------------------------------------------------------
     * SOLVE
     * --------------------------------------------------------
     */
    harmonic_solution_init(
        solution
    );

    if (harmonic_solve(
        &harmonic_system,
        &model->P,
        &excitation,
        solution,
        config->parallel_threads
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: harmonic PEEC solve failed.\n"
        );

        harmonic_solution_free(
            solution
        );

        harmonic_system_free(
            &harmonic_system
        );

        free(
            U_inc
        );

        return -1;
    }


    /*
     * --------------------------------------------------------
     * RESIDUAL
     * --------------------------------------------------------
     */
    HarmonicResidual residual;

    if (harmonic_compute_residual(
        &harmonic_system,
        &model->P,
        &excitation,
        solution,
        config->parallel_threads,
        &residual
    ) != 0)
    {
        harmonic_solution_free(
            solution
        );

        harmonic_system_free(
            &harmonic_system
        );

        free(
            U_inc
        );

        return -1;
    }

    printf("\n");
    printf("============================================================\n");
    printf("HARMONIC PEEC SOLUTION\n");
    printf("============================================================\n");

    printf(
        "Frequency          : %.9e Hz\n",
        frequency
    );

    printf(
        "Maximum |I_e|      : %.9e A\n",
        complex_vector_max_abs(
            solution->edge_current,
            solution->n_edges
        )
    );

    printf(
        "Maximum |V_j|      : %.9e V\n",
        complex_vector_max_abs(
            solution->node_voltage,
            solution->n_nodes
        )
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

    printf("============================================================\n");

    int result =
    0;

    if (!isfinite(
        residual.backward_error
    ) ||
    residual.backward_error > 1.0e-10)
    {
        fprintf(
            stderr,
            "ERROR: harmonic backward error is too large.\n"
        );

        result =
        -1;
    }

    harmonic_system_free(
        &harmonic_system
    );

    free(
        U_inc
    );

    if (result != 0)
    {
        harmonic_solution_free(
            solution
        );
    }

    return result;
}


/*
 * ============================================================
 * HARMONIC PARAVIEW
 * ============================================================
 */
static int write_harmonic_vtk_if_requested(
    const PeecConfig *config,
    const PeecModel *model,
    const char *prefix,
    const HarmonicSolution *solution
)
{
    if (config == NULL ||
        model == NULL ||
        prefix == NULL ||
        solution == NULL)
    {
        return -1;
    }

    if (!config->write_vtk)
    {
        return 0;
    }

    ChargeRecovery recovery;

    charge_recovery_init(
        &recovery
    );

    if (charge_recovery_factorize(
        &recovery,
        &model->P,
        config->parallel_threads
    ) != 0)
    {
        charge_recovery_free(
            &recovery
        );

        return -1;
    }

    Complex *Q =
    calloc(
        model->mesh.n_nodes,
        sizeof(Complex)
    );

    if (Q == NULL)
    {
        charge_recovery_free(
            &recovery
        );

        return -1;
    }

    if (charge_recovery_complex(
        &recovery,
        solution->node_voltage,
        Q
    ) != 0)
    {
        free(
            Q
        );

        charge_recovery_free(
            &recovery
        );

        return -1;
    }

    int result =
    paraview_write_complex_state(
        config->vtk_directory,
        prefix,
        &model->mesh,
        &model->dual,
        solution->node_voltage,
        Q,
        solution->edge_current
    );

    free(
        Q
    );

    charge_recovery_free(
        &recovery
    );

    return result;
}


/*
 * ============================================================
 * RCS CSV
 * ============================================================
 */
static int save_rcs_sweep(
    const PeecConfig *config,
    double frequency,
    double wave_number,
    const PeecModel *model,
    const IncidentField *field,
    const HarmonicSolution *solution
)
{
    if (config == NULL ||
        model == NULL ||
        field == NULL ||
        solution == NULL)
    {
        return -1;
    }

    if (ensure_parent_directory(
        config->output_file
    ) != 0)
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
            "ERROR: cannot open RCS output '%s'.\n",
            config->output_file
        );

        return -1;
    }

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

    RcsOptions options;

    rcs_options_default(
        &options,
        config->parallel_threads
    );

    options.integration =
    config->rcs_use_dual
    ? RCS_CURRENT_DUAL_PATCH_QUADRATURE
    : RCS_CURRENT_MIDPOINT;

    options.surface_order =
    config->rcs_order;


    /*
     * --------------------------------------------------------
     * BACKSCATTER
     * --------------------------------------------------------
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

    RcsResult backscatter;

    if (rcs_compute_direction(
        &model->mesh,
        &model->dual,
        field,
        solution,
        &options,
        back_theta,
        back_phi,
        &backscatter
    ) != 0)
    {
        fclose(
            file
        );

        return -1;
    }


    /*
     * --------------------------------------------------------
     * SWEEP
     * --------------------------------------------------------
     */
    double max_sigma =
    -1.0;

    double max_dbsm =
    -INFINITY;

    double max_phi =
    0.0;

    size_t count =
    0;

    for (double phi = config->rcs_phi_start;
         phi <= config->rcs_phi_end + 0.5 * config->rcs_phi_step;
    phi += config->rcs_phi_step)
         {
             RcsResult result;

             if (rcs_compute_direction(
                 &model->mesh,
                 &model->dual,
                 field,
                 solution,
                 &options,
                 config->rcs_observation_theta,
                 phi,
                 &result
             ) != 0)
             {
                 fclose(
                     file
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
                 frequency,
                 wave_number,
                 config->ka,
                 config->characteristic_length
             );

             if (result.sigma >
                 max_sigma)
             {
                 max_sigma =
                 result.sigma;

                 max_dbsm =
                 result.sigma_dbsm;

                 max_phi =
                 phi;
             }

             ++count;
         }

         fclose(
             file
         );

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
             backscatter.sigma
         );

         printf(
             "Backscatter RCS    : %.9f dBsm\n",
             backscatter.sigma_dbsm
         );

         printf(
             "Maximum sweep RCS  : %.9e m^2\n",
             max_sigma
         );

         printf(
             "Maximum sweep RCS  : %.9f dBsm\n",
             max_dbsm
         );

         printf(
             "Maximum at phi     : %.6f deg\n",
             max_phi
         );

         printf(
             "Sweep points       : %zu\n",
             count
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
 * LIGHTNING SOURCE CONFIGURATION
 * ============================================================
 */
static int configure_lightning_source(
    const PeecConfig *config,
    const Mesh *mesh,
    LightningSource *source
)
{
    if (config == NULL ||
        mesh == NULL ||
        source == NULL)
    {
        return -1;
    }

    lightning_source_default(
        source
    );

    source->K =
    config->lightning_K;

    source->alpha =
    config->lightning_alpha;

    source->beta =
    config->lightning_beta;

    source->delay =
    config->lightning_delay;

    double strike_error =
    0.0;

    double return_error =
    0.0;

    if (config->strike_node_set &&
        config->return_node_set)
    {
        if ((size_t)config->strike_node >= mesh->n_nodes ||
            (size_t)config->return_node >= mesh->n_nodes)
        {
            fprintf(
                stderr,
                "ERROR: lightning node index is outside mesh: "
                "strike=%d, return=%d, n_nodes=%zu.\n",
                config->strike_node,
                config->return_node,
                mesh->n_nodes
            );

            return -1;
        }

        source->strike_node =
        (size_t)config->strike_node;

        source->return_node =
        (size_t)config->return_node;
    }
    else
    {
        if (lightning_find_nearest_node(
            mesh,
            config->strike_xyz,
            &source->strike_node,
            &strike_error
        ) != 0)
        {
            return -1;
        }

        if (lightning_find_nearest_node(
            mesh,
            config->return_xyz,
            &source->return_node,
            &return_error
        ) != 0)
        {
            return -1;
        }
    }

    if (lightning_source_validate(
        source,
        mesh->n_nodes
    ) != 0)
    {
        return -1;
    }

    printf("\n");
    printf("Lightning node mapping\n");
    printf("------------------------------------------------------------\n");

    if (config->strike_node_set &&
        config->return_node_set)
    {
        printf(
            "Strike node        : %zu (direct index)\n",
               source->strike_node
        );

        printf(
            "Return node        : %zu (direct index)\n",
               source->return_node
        );
    }
    else
    {
        printf(
            "Strike node        : %zu, error %.9e m\n",
            source->strike_node,
            strike_error
        );

        printf(
            "Return node        : %zu, error %.9e m\n",
            source->return_node,
            return_error
        );
    }

    lightning_print_info(
        source,
        mesh
    );

    return 0;
}


/*
 * ============================================================
 * NODE KCL DIAGNOSTICS
 * ============================================================
 *
 * Для новой transient-формы ЧЭС:
 *
 *     A^T I - F dV_c/dt + J = 0,
 *
 * где:
 *
 *     F_j = 1/P_jj.
 *
 * Функция ниже вычисляет:
 *
 *     (A^T I)_j
 *
 * непосредственно из матрицы инцидентности.
 */
static double node_incidence_current(
    const IncidenceMatrix *A,
    const double *edge_current,
    size_t node
)
{
    if (A == NULL ||
        A->data == NULL ||
        edge_current == NULL ||
        node >= A->cols)
    {
        return NAN;
    }

    double value =
    0.0;

    for (size_t e = 0;
         e < A->rows;
    ++e)
         {
             value +=
             A->data[
                 e * A->cols
                 + node
             ]
             * edge_current[e];
         }

         return value;
}


static double edge_current_max_abs(
    const double *edge_current,
    size_t n_edges,
    size_t *edge_index
)
{
    if (edge_current == NULL ||
        n_edges == 0)
    {
        if (edge_index != NULL)
        {
            *edge_index =
            0;
        }

        return 0.0;
    }

    size_t index =
    0;

    double maximum =
    fabs(
        edge_current[0]
    );

    for (size_t e = 1;
         e < n_edges;
    ++e)
         {
             double value =
             fabs(
                 edge_current[e]
             );

             if (value > maximum)
             {
                 maximum =
                 value;

                 index =
                 e;
             }
         }

         if (edge_index != NULL)
         {
             *edge_index =
             index;
         }

         return maximum;
}


/*
 * ============================================================
 * STAGE-2 SHUNT VISUALIZATION METADATA
 * ============================================================
 *
 * Для каждого сохраненного VTK-кадра Stage 2 создается небольшой
 * sidecar-файл:
 *
 *     <prefix>_shunt_XXXXXX.dat
 *
 * Он содержит:
 *
 *     time
 *     node_a, xyz_a
 *     node_b, xyz_b
 *     I_shunt
 *     R_shunt
 *
 * Это позволяет визуализатору автоматически показать точки A/B
 * и текущее значение тока шунта, не меняя формат основного VTK.
 */
static int write_shunt_visualization_metadata(
    const PeecConfig *config,
    const char *prefix,
    size_t frame,
    double time,
    const Mesh *mesh,
    const InternalCircuitShunt *shunt,
    double shunt_current
)
{
    if (config == NULL ||
        prefix == NULL ||
        mesh == NULL ||
        shunt == NULL ||
        !config->write_vtk)
    {
        return -1;
    }

    if (shunt->node_a >= mesh->n_nodes ||
        shunt->node_b >= mesh->n_nodes)
    {
        return -1;
    }

    if (ensure_directory(
        config->vtk_directory
    ) != 0)
    {
        return -1;
    }

    char filename[4096];

    int written =
    snprintf(
        filename,
        sizeof(filename),
             "%s/%s_shunt_%06zu.dat",
             config->vtk_directory,
             prefix,
             frame
    );

    if (written < 0 ||
        (size_t)written >= sizeof(filename))
    {
        return -1;
    }

    FILE *file =
    fopen(
        filename,
        "w"
    );

    if (file == NULL)
    {
        return -1;
    }

    const double *a =
    &mesh->xyz[
        3 * shunt->node_a
    ];

    const double *b =
    &mesh->xyz[
        3 * shunt->node_b
    ];

    fprintf(
        file,
        "time %.17e\n",
        time
    );

    fprintf(
        file,
        "node_a %zu\n",
        shunt->node_a
    );

    fprintf(
        file,
        "xyz_a %.17e %.17e %.17e\n",
        a[0],
        a[1],
        a[2]
    );

    fprintf(
        file,
        "node_b %zu\n",
        shunt->node_b
    );

    fprintf(
        file,
        "xyz_b %.17e %.17e %.17e\n",
        b[0],
        b[1],
        b[2]
    );

    fprintf(
        file,
        "shunt_current_A %.17e\n",
        shunt_current
    );

    fprintf(
        file,
        "shunt_resistance_ohm %.17e\n",
        shunt->resistance
    );

    fclose(
        file
    );

    return 0;
}


/*
 * ============================================================
 * SLOT-CELL VISUALIZATION METADATA
 * ============================================================
 *
 * Для каждого сохраненного transient-кадра создается:
 *
 *     <prefix>_slot_XXXXXX.dat
 *
 * Файл содержит ток индуктивной части каждой эквивалентной
 * ячейки щели, а также ее параметры L, C и длину поддержки.
 */
static int write_slot_visualization_metadata(
    const PeecConfig *config,
    const char *prefix,
    size_t frame,
    double time,
    const SlotCircuit *slot,
    const double *slot_current
)
{
    if (config == NULL ||
        prefix == NULL ||
        slot == NULL ||
        slot_current == NULL ||
        !config->write_vtk)
    {
        return -1;
    }

    if (ensure_directory(
        config->vtk_directory
    ) != 0)
    {
        return -1;
    }

    char filename[4096];

    int written =
    snprintf(
        filename,
        sizeof(filename),
             "%s/%s_slot_%06zu.dat",
             config->vtk_directory,
             prefix,
             frame
    );

    if (written < 0 ||
        (size_t)written >= sizeof(filename))
    {
        return -1;
    }

    FILE *file =
    fopen(
        filename,
        "w"
    );

    if (file == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot open slot metadata file '%s'.\n",
            filename
        );

        return -1;
    }

    fprintf(
        file,
        "time %.17e\n",
        time
    );

    fprintf(
        file,
        "n_cells %zu\n",
        slot->n_cells
    );

    for (size_t m = 0;
         m < slot->n_cells;
    ++m)
         {
             const SlotCell *cell =
             &slot->cells[m];

             fprintf(
                 file,
                 "cell %zu "
                 "node_a %zu "
                 "node_b %zu "
                 "I_L_A %.17e "
                 "L_H %.17e "
                 "C_F %.17e "
                 "dl_m %.17e\n",
                 m,
                 cell->node_a,
                 cell->node_b,
                 slot_current[m],
                 cell->inductance,
                 cell->capacitance,
                 cell->support_length
             );
         }

         fclose(
             file
         );

         return 0;
}


/*
 * ============================================================
 * APERTURE CURRENT TRANSFER DIAGNOSTIC
 * ============================================================
 *
 * Prints the exact sparse transfer
 *
 *     I_gamma = B_gamma I_cover
 *
 * at one stored time sample.
 */
static void print_aperture_current_transfer_diagnostic(
    const ApertureCoupling *coupling,
    const Mesh *closed_mesh,
    const double *closed_edge_current,
    const double *open_node_current,
    double time
)
{
    if (coupling == NULL ||
        closed_mesh == NULL ||
        closed_mesh->edges == NULL ||
        closed_mesh->xyz == NULL ||
        closed_edge_current == NULL ||
        open_node_current == NULL)
    {
        return;
    }

    printf("\n");
    printf("============================================================\n");
    printf("APERTURE CURRENT TRANSFER DIAGNOSTIC\n");
    printf("============================================================\n");
    printf(
        "time                : %.9e s (%.6f us)\n",
           time,
           time * 1.0e6
    );

    printf(
        "cut edges           : %zu\n",
        coupling->n_cut_edges
    );

    printf("\n");
    printf(
        " k   edge   boundary  interior  open   sign"
        "        I_edge [A]        contribution [A]\n"
    );

    printf(
        "------------------------------------------------------------"
        "----------------------\n"
    );

    for (size_t k = 0;
         k < coupling->n_cut_edges;
    ++k)
         {
             const ApertureCutEdge *cut =
             &coupling->cut_edges[k];

             size_t edge =
             cut->closed_edge;

             size_t a =
             (size_t)closed_mesh->edges[
                 2 * edge + 0
             ];

             size_t b =
             (size_t)closed_mesh->edges[
                 2 * edge + 1
             ];

             size_t boundary_node =
             cut->sign < 0.0
             ? a
             : b;

             size_t interior_node =
             cut->sign < 0.0
             ? b
             : a;

             double I_edge =
             closed_edge_current[edge];

             double contribution =
             cut->sign
             * I_edge;

             printf(
                 "%2zu  %5zu  %8zu  %8zu  %4zu  %+4.0f"
                 "  % .12e  % .12e\n",
                 k,
                 edge,
                 boundary_node,
                 interior_node,
                 cut->open_node,
                 cut->sign,
                 I_edge,
                 contribution
             );
         }

         printf("\n");
         printf("OPEN-NODE I_gamma AFTER AGGREGATION\n");
         printf("------------------------------------------------------------\n");
         printf(
             " open node        I_gamma [A]\n"
         );

         double net =
         0.0;

         double half_l1 =
         0.0;

         double max_abs =
         0.0;

         size_t max_node =
         0;

         size_t nonzero =
         0;

         for (size_t j = 0;
              j < coupling->n_open_nodes;
    ++j)
              {
                  double value =
                  open_node_current[j];

                  net +=
                  value;

                  half_l1 +=
                  fabs(
                      value
                  );

                  if (fabs(value) >
                      max_abs)
                  {
                      max_abs =
                      fabs(
                          value
                      );

                      max_node =
                      j;
                  }

                  /*
                   * Print only actually excited nodes.
                   */
                  if (fabs(value) >
                      1.0e-12)
                  {
                      printf(
                          " %8zu   % .12e\n",
                          j,
                          value
                      );

                      ++nonzero;
                  }
              }

              half_l1 *=
              0.5;

              printf("\n");
              printf(
                  "nonzero I_gamma nodes : %zu\n",
                  nonzero
              );

              printf(
                  "net I_gamma           : %.12e A\n",
                  net
              );

              printf(
                  "max |I_gamma|         : %.12e A at open node %zu\n",
                  max_abs,
                  max_node
              );

              printf(
                  "0.5 sum |I_gamma|     : %.12e A\n",
                  half_l1
              );

              printf("============================================================\n");
}


/*
 * ============================================================
 * TRANSIENT VTK OUTPUT
 * ============================================================
 */
static int write_transient_vtk(
    const PeecConfig *config,
    const char *prefix,
    size_t frame,
    double time,
    const PeecModel *model,
    const TransientSystem *transient,
    const double *I,
    const double *V_c,
    double *phi,
    double *Q
)
{
    if (config == NULL ||
        prefix == NULL ||
        model == NULL ||
        transient == NULL ||
        I == NULL ||
        V_c == NULL)
    {
        return -1;
    }

    if (!config->write_vtk)
    {
        return 0;
    }

    if (phi == NULL ||
        Q == NULL)
    {
        return -1;
    }

    if (transient_recover_node_state(
        transient,
        V_c,
        phi,
        Q
    ) != 0)
    {
        return -1;
    }

    /*
     * VTK получает физический полный потенциал:
     *
     *     phi = S V_c,
     *
     * а не внутреннюю неизвестную V_c.
     */
    return paraview_write_real_state(
        config->vtk_directory,
        prefix,
        frame,
        time,
        &model->mesh,
        &model->dual,
        phi,
        Q,
        I
    );
}


/*
 * ============================================================
 * MESH-INFO
 * ============================================================
 */
static int run_mesh_info(
    const PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    PeecModel model;

    peec_model_init(
        &model
    );

    int result =
    peec_load_mesh_geometry(
        config->mesh_file,
        &model
    );

    if (result == 0)
    {
        long long euler =
        (long long)model.mesh.n_nodes
        - (long long)model.mesh.n_edges
        + (long long)model.mesh.n_quads;

        printf(
            "Euler characteristic : %lld\n",
            euler
        );
    }

    peec_model_free(
        &model
    );

    return result;
}


/*
 * ============================================================
 * SCATTERING
 * ============================================================
 */
static int run_scattering(
    const PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    PeecModel model;

    peec_model_init(
        &model
    );

    HarmonicSolution solution;

    harmonic_solution_init(
        &solution
    );

    IncidentField field;

    incident_field_init(
        &field
    );

    int result =
    -1;

    if (peec_assemble_model(
        config->mesh_file,
        config->parallel_threads,
        &model
    ) != 0)
    {
        goto cleanup;
    }

    if (peec_export_matrices_if_requested(
        config,
        &model,
        config->matrix_directory
    ) != 0)
    {
        goto cleanup;
    }

    if (solve_harmonic_scattering(
        config,
        &model,
        config->frequency,
        &field,
        &solution
    ) != 0)
    {
        goto cleanup;
    }

    if (write_harmonic_vtk_if_requested(
        config,
        &model,
        "scattering",
        &solution
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: ParaView scattering output failed.\n"
        );

        goto cleanup;
    }

    result =
    0;


    cleanup:

    harmonic_solution_free(
        &solution
    );

    peec_model_free(
        &model
    );

    return result;
}


/*
 * ============================================================
 * RCS
 * ============================================================
 */
static int run_rcs(
    const PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    /*
     * Для RCS основной способ задания частоты:
     *
     *     ka = k a.
     */
    double wave_number =
    config->ka
    / config->characteristic_length;

    double frequency =
    wave_number
    * C0
    / (
        2.0
        * M_PI
    );

    if (!isfinite(
        frequency
    ) ||
    frequency <= 0.0)
    {
        return -1;
    }

    printf("\n");
    printf("============================================================\n");
    printf("RCS ELECTRICAL SIZE\n");
    printf("============================================================\n");

    printf(
        "a                  : %.9e m\n",
        config->characteristic_length
    );

    printf(
        "ka                 : %.9e\n",
        config->ka
    );

    printf(
        "k                  : %.9e rad/m\n",
        wave_number
    );

    printf(
        "Frequency          : %.9e Hz\n",
        frequency
    );

    printf("============================================================\n");


    PeecModel model;

    peec_model_init(
        &model
    );

    HarmonicSolution solution;

    harmonic_solution_init(
        &solution
    );

    IncidentField field;

    incident_field_init(
        &field
    );

    int result =
    -1;

    if (peec_assemble_model(
        config->mesh_file,
        config->parallel_threads,
        &model
    ) != 0)
    {
        goto cleanup;
    }

    if (peec_export_matrices_if_requested(
        config,
        &model,
        config->matrix_directory
    ) != 0)
    {
        goto cleanup;
    }

    if (solve_harmonic_scattering(
        config,
        &model,
        frequency,
        &field,
        &solution
    ) != 0)
    {
        goto cleanup;
    }

    if (save_rcs_sweep(
        config,
        frequency,
        wave_number,
        &model,
        &field,
        &solution
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: RCS sweep failed.\n"
        );

        goto cleanup;
    }

    if (write_harmonic_vtk_if_requested(
        config,
        &model,
        "rcs_solution",
        &solution
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: ParaView RCS solution output failed.\n"
        );

        goto cleanup;
    }

    result =
    0;


    cleanup:

    harmonic_solution_free(
        &solution
    );

    peec_model_free(
        &model
    );

    return result;
}


/*
 * ============================================================
 * LIGHTNING CLOSED BODY
 * ============================================================
 */
static int run_lightning_closed(
    const PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    PeecModel model;

    peec_model_init(
        &model
    );

    TransientSystem transient;

    transient_system_init(
        &transient
    );

    LightningSource source;

    IncidentField incident_pulse;

    incident_field_init(
        &incident_pulse
    );

    ScatteringOptions scattering_options;

    scattering_options_default(
        &scattering_options,
        config->parallel_threads
    );

    double *I =
    NULL;

    double *V =
    NULL;

    double *phi =
    NULL;

    double *Q =
    NULL;

    double *source_n =
    NULL;

    double *source_np1 =
    NULL;

    double *edge_voltage_n =
    NULL;

    double *edge_voltage_np1 =
    NULL;

    int result =
    -1;

    if (peec_assemble_model(
        config->closed_mesh_file,
        config->parallel_threads,
        &model
    ) != 0)
    {
        goto cleanup;
    }

    if (peec_export_matrices_if_requested(
        config,
        &model,
        config->matrix_directory
    ) != 0)
    {
        goto cleanup;
    }

    if (config->excitation ==
        PEEC_EXCITATION_LIGHTNING_CURRENT)
    {
        if (configure_lightning_source(
            config,
            &model.mesh,
            &source
        ) != 0)
        {
            fprintf(
                stderr,
                "ERROR: lightning source configuration failed.\n"
            );

            goto cleanup;
        }
    }
    else
    {
        if (configure_incident_pulse(
            config,
            &incident_pulse
        ) != 0)
        {
            goto cleanup;
        }
    }

    if (transient_system_create_stage1(
        &transient,
        &model.L,
        &model.R,
        &model.incidence,
        &model.P,
        transient_scheme_from_config(
            config
        ),
        config->dt
    ) != 0)
    {
        goto cleanup;
    }

    if (transient_factorize(
        &transient,
        config->parallel_threads
    ) != 0)
    {
        goto cleanup;
    }

    transient_print_info(
        &transient
    );

    I = calloc(
        model.mesh.n_edges,
        sizeof(double)
    );

    V = calloc(
        model.mesh.n_nodes,
        sizeof(double)
    );

    source_n = calloc(
        model.mesh.n_nodes,
        sizeof(double)
    );

    source_np1 = calloc(
        model.mesh.n_nodes,
        sizeof(double)
    );

    edge_voltage_n = calloc(
        model.mesh.n_edges,
        sizeof(double)
    );

    edge_voltage_np1 = calloc(
        model.mesh.n_edges,
        sizeof(double)
    );

    if (I == NULL ||
        V == NULL ||
        source_n == NULL ||
        source_np1 == NULL ||
        edge_voltage_n == NULL ||
        edge_voltage_np1 == NULL)
    {
        goto cleanup;
    }

    if (config->write_vtk)
    {
        phi = calloc(
            model.mesh.n_nodes,
            sizeof(double)
        );

        Q = calloc(
            model.mesh.n_nodes,
            sizeof(double)
        );

        if (phi == NULL ||
            Q == NULL)
        {
            goto cleanup;
        }

        if (write_transient_vtk(
            config,
            "lightning_closed",
            0,
            0.0,
            &model,
            &transient,
            I,
            V,
            phi,
            Q
        ) != 0)
        {
            goto cleanup;
        }
    }

    size_t n_steps =
    0;

    if (compute_constant_dt_steps(
        config->t_end,
        config->dt,
        &n_steps
    ) != 0)
    {
        goto cleanup;
    }

    size_t frame =
    1;

    size_t progress_step =
    n_steps / 20;

    if (progress_step < 1)
    {
        progress_step =
        1;
    }

    for (size_t step = 0;
         step < n_steps;
    ++step)
         {
             double t_n =
             step
             * config->dt;

             double t_np1 =
             (step + 1)
             * config->dt;

             if (t_np1 >
                 config->t_end)
             {
                 t_np1 =
                 config->t_end;
             }

             /*
              * Последний укороченный шаг не поддерживаем:
              *
              * transient matrix построена для постоянного dt.
              *
              * Поэтому t_end лучше выбирать кратным dt.
              */
             if (fabs(
                 t_np1 - t_n - config->dt
             ) >
             1.0e-12
             * fmax(
                 1.0,
                 config->dt
             ))
             {
                 fprintf(
                     stderr,
                     "ERROR: t_end must be an integer multiple of dt "
                     "for the current constant-dt transient solver.\n"
                 );

                 goto cleanup;
             }

             const double *node_source_n_ptr =
             NULL;

             const double *node_source_np1_ptr =
             NULL;

             const double *edge_voltage_n_ptr =
             NULL;

             const double *edge_voltage_np1_ptr =
             NULL;

             if (config->excitation ==
                 PEEC_EXCITATION_LIGHTNING_CURRENT)
             {
                 if (lightning_build_node_current(
                     &source,
                     model.mesh.n_nodes,
                     t_n,
                     source_n
                 ) != 0 ||
                 lightning_build_node_current(
                     &source,
                     model.mesh.n_nodes,
                     t_np1,
                     source_np1
                 ) != 0)
                 {
                     goto cleanup;
                 }

                 node_source_n_ptr =
                 source_n;

                 node_source_np1_ptr =
                 source_np1;
             }
             else
             {
                 if (scattering_build_time_excitation(
                     &model.mesh,
                     &incident_pulse,
                     &scattering_options,
                     t_n,
                     edge_voltage_n
                 ) != 0 ||
                 scattering_build_time_excitation(
                     &model.mesh,
                     &incident_pulse,
                     &scattering_options,
                     t_np1,
                     edge_voltage_np1
                 ) != 0)
                 {
                     fprintf(
                         stderr,
                         "ERROR: cannot build incident-pulse edge excitation.\n"
                     );

                     goto cleanup;
                 }

                 edge_voltage_n_ptr =
                 edge_voltage_n;

                 edge_voltage_np1_ptr =
                 edge_voltage_np1;
             }

             if (transient_build_rhs_stage1(
                 &transient,
                 I,
                 V,
                 node_source_n_ptr,
                 node_source_np1_ptr,
                 edge_voltage_n_ptr,
                 edge_voltage_np1_ptr
             ) != 0 ||
             transient_solve_step(
                 &transient
             ) != 0 ||
             transient_extract_stage1_solution(
                 &transient,
                 I,
                 V
             ) != 0)
             {
                 goto cleanup;
             }

             if (config->write_vtk &&
                 (
                     (step + 1) % (size_t)config->vtk_every == 0 ||
                     step + 1 == n_steps
                 ))
             {
                 if (write_transient_vtk(
                     config,
                     "lightning_closed",
                     frame++,
                     t_np1,
                     &model,
                     &transient,
                     I,
                     V,
                     phi,
                     Q
                 ) != 0)
                 {
                     goto cleanup;
                 }
             }

             if ((step + 1) % progress_step == 0 ||
                 step + 1 == n_steps)
             {
                 if (config->excitation ==
                     PEEC_EXCITATION_LIGHTNING_CURRENT)
                 {
                     printf(
                         "Transient: %6.1f%%  t = %.9e s  I_L = %.9e A\n",
                         100.0
                         * (double)(step + 1)
                         / (double)n_steps,
                            t_np1,
                            lightning_current(
                                &source,
                                t_np1
                            )
                     );
                 }
                 else
                 {
                     printf(
                         "Transient: %6.1f%%  t = %.9e s  incident-pulse\n",
                         100.0
                         * (double)(step + 1)
                         / (double)n_steps,
                            t_np1
                     );
                 }
             }
         }

         result =
         0;


         cleanup:

         free(
             edge_voltage_np1
         );

         free(
             edge_voltage_n
         );

         free(
             source_np1
         );

         free(
             source_n
         );

         free(
             Q
         );

         free(
             phi
         );

         free(
             V
         );

         free(
             I
         );

         transient_system_free(
             &transient
         );

         peec_model_free(
             &model
         );

         return result;
}


/*
 * ============================================================
 * LIGHTNING TWO-STAGE
 * ============================================================
 *
 * Stage 1 и Stage 2 интегрируются синхронно по одной
 * временной сетке.
 *
 *
 * Это позволяет НЕ хранить всю историю:
 *
 *     I_cover(t).
 *
 *
 * На каждом шаге:
 *
 *     1. решаем Stage 1;
 *
 *     2. строим
 *
 *            I_gamma^n
 *            I_gamma^(n+1);
 *
 *     3. сразу решаем Stage 2.
 *
 *
 * Математически это эквивалентно последовательному запуску,
 * поскольку Stage 1 не зависит от Stage 2.
 */
static int run_lightning_two_stage(
    const PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    PeecModel closed_model;
    PeecModel open_model;

    peec_model_init(
        &closed_model
    );

    peec_model_init(
        &open_model
    );

    ApertureCoupling coupling;

    aperture_coupling_init(
        &coupling
    );

    LightningSource lightning;

    IncidentField incident_pulse;

    incident_field_init(
        &incident_pulse
    );

    ScatteringOptions scattering_options;

    scattering_options_default(
        &scattering_options,
        config->parallel_threads
    );

    InternalCircuitShunt shunt;

    internal_circuit_shunt_default(
        &shunt
    );

    TransientSystem stage1;
    TransientSystem stage2;

    transient_system_init(
        &stage1
    );

    transient_system_init(
        &stage2
    );

    double *I1 =
    NULL;

    double *V1 =
    NULL;

    double *phi1 =
    NULL;

    double *Q1 =
    NULL;

    double *lightning_n =
    NULL;

    double *lightning_np1 =
    NULL;

    double *U1_n =
    NULL;

    double *U1_np1 =
    NULL;

    double *I2 =
    NULL;

    double *V2 =
    NULL;

    double *phi2 =
    NULL;

    double *Q2 =
    NULL;

    double *I_gamma_n =
    NULL;

    double *I_gamma_np1 =
    NULL;

    double I_shunt =
    0.0;

    int stage2_has_shunt =
    config->use_internal_shunt;

    int result =
    -1;

    /*
     * --------------------------------------------------------
     * STAGE-1 KCL / MAXIMUM CURRENT DIAGNOSTICS
     * --------------------------------------------------------
     */
    double sampled_source_max =
    -1.0;

    double sampled_source_time =
    0.0;

    double sampled_strike_ATI =
    0.0;

    double sampled_strike_Cdot =
    0.0;

    double sampled_strike_J =
    0.0;

    double sampled_strike_residual =
    0.0;

    double sampled_return_ATI =
    0.0;

    double sampled_return_Cdot =
    0.0;

    double sampled_return_J =
    0.0;

    double sampled_return_residual =
    0.0;

    double stage1_edge_max =
    0.0;

    double stage1_edge_max_time =
    0.0;

    size_t stage1_edge_max_index =
    0;

    /*
     * --------------------------------------------------------
     * STAGE-2 APERTURE / PORT DIAGNOSTICS
     * --------------------------------------------------------
     *
     * These diagnostics do not alter the equations.
     *
     * We track:
     *
     *   max_j |I_gamma,j|
     *   0.5 * sum_j |I_gamma,j|
     *   V_AB = phi_A - phi_B
     *   max |I_shunt| when the shunt is enabled.
     *
     * The A/B nodes may also be supplied without --shunt-R.
     * In that case Stage 2 is open-circuit and V_AB is V_oc.
     */
    int stage2_probe_enabled =
    0;

    size_t stage2_probe_a =
    0;

    size_t stage2_probe_b =
    0;

    double gamma_max_global =
    0.0;

    double gamma_max_global_time =
    0.0;

    size_t gamma_max_global_node =
    0;

    double gamma_half_l1_at_global_max =
    0.0;

    double gamma_half_l1_global =
    0.0;

    double gamma_half_l1_global_time =
    0.0;

    double port_voltage_abs_max =
    0.0;

    double port_voltage_at_abs_max =
    0.0;

    double port_voltage_abs_max_time =
    0.0;

    double shunt_current_abs_max =
    0.0;

    double shunt_current_at_abs_max =
    0.0;

    double shunt_current_abs_max_time =
    0.0;

    /*
     * Snapshot of aperture transfer at the sampled source peak.
     */
    double *aperture_peak_edge_current =
    NULL;

    double *aperture_peak_gamma =
    NULL;

    double aperture_peak_time =
    0.0;


    /*
     * --------------------------------------------------------
     * BOTH PEEC MODELS
     * --------------------------------------------------------
     */
    printf("\n");
    printf("============================================================\n");
    printf("LIGHTNING STAGE 1 / CLOSED MODEL\n");
    printf("============================================================\n");

    if (peec_assemble_model(
        config->closed_mesh_file,
        config->parallel_threads,
        &closed_model
    ) != 0)
    {
        goto cleanup;
    }

    printf("\n");
    printf("============================================================\n");
    printf("LIGHTNING STAGE 2 / OPEN MODEL\n");
    printf("============================================================\n");

    if (peec_assemble_model(
        config->open_mesh_file,
        config->parallel_threads,
        &open_model
    ) != 0)
    {
        goto cleanup;
    }


    /*
     * --------------------------------------------------------
     * MATRIX EXPORT
     * --------------------------------------------------------
     */
    if (config->save_matrices)
    {
        char stage1_dir[4096];
        char stage2_dir[4096];

        if (make_subdirectory(
            stage1_dir,
            sizeof(stage1_dir),
                              config->matrix_directory,
                              "stage1"
        ) != 0 ||
        make_subdirectory(
            stage2_dir,
            sizeof(stage2_dir),
                          config->matrix_directory,
                          "stage2"
        ) != 0)
        {
            goto cleanup;
        }

        if (matrix_export_peec_set(
            stage1_dir,
            &closed_model.mesh,
            &closed_model.L,
            &closed_model.P,
            &closed_model.PAT,
            &closed_model.R
        ) != 0 ||
        matrix_export_peec_set(
            stage2_dir,
            &open_model.mesh,
            &open_model.L,
            &open_model.P,
            &open_model.PAT,
            &open_model.R
        ) != 0)
        {
            goto cleanup;
        }
    }


    /*
     * --------------------------------------------------------
     * APERTURE COUPLING
     * --------------------------------------------------------
     */
    if (aperture_coupling_load_json(
        &coupling,
        &closed_model.mesh,
        &open_model.mesh,
        config->aperture_map_file
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: aperture coupling setup failed.\n"
        );

        goto cleanup;
    }

    aperture_coupling_print_info(
        &coupling
    );


    /*
     * --------------------------------------------------------
     * LIGHTNING SOURCE
     * --------------------------------------------------------
     */
    if (config->excitation ==
        PEEC_EXCITATION_LIGHTNING_CURRENT)
    {
        if (configure_lightning_source(
            config,
            &closed_model.mesh,
            &lightning
        ) != 0)
        {
            goto cleanup;
        }
    }
    else
    {
        if (configure_incident_pulse(
            config,
            &incident_pulse
        ) != 0)
        {
            goto cleanup;
        }
    }


    /*
     * --------------------------------------------------------
     * STAGE 1 TRANSIENT
     * --------------------------------------------------------
     */
    if (transient_system_create_stage1(
        &stage1,
        &closed_model.L,
        &closed_model.R,
        &closed_model.incidence,
        &closed_model.P,
        transient_scheme_from_config(
            config
        ),
        config->dt
    ) != 0 ||
    transient_factorize(
        &stage1,
        config->parallel_threads
    ) != 0)
    {
        goto cleanup;
    }


    /*
     * --------------------------------------------------------
     * STAGE 2 TRANSIENT
     * --------------------------------------------------------
     *
     * Если внутреннего шунта нет, уравнения Stage 2 по форме
     * совпадают с обычной Stage-1 системой, только источником
     * является:
     *
     *     I_gamma.
     */
    if (stage2_has_shunt)
    {
        shunt.resistance =
        config->shunt_resistance;

        if (config->shunt_a_node_set &&
            config->shunt_b_node_set)
        {
            if ((size_t)config->shunt_a_node >= open_model.mesh.n_nodes ||
                (size_t)config->shunt_b_node >= open_model.mesh.n_nodes)
            {
                fprintf(
                    stderr,
                    "ERROR: shunt node index is outside Stage-2 mesh: "
                    "A=%d, B=%d, n_nodes=%zu.\n",
                    config->shunt_a_node,
                    config->shunt_b_node,
                    open_model.mesh.n_nodes
                );

                goto cleanup;
            }

            if (internal_circuit_shunt_set_nodes(
                &shunt,
                &open_model.mesh,
                (size_t)config->shunt_a_node,
                                                 (size_t)config->shunt_b_node
            ) != 0)
            {
                goto cleanup;
            }
        }
        else
        {
            if (internal_circuit_shunt_map_coordinates(
                &shunt,
                &open_model.mesh,
                config->shunt_a_xyz,
                config->shunt_b_xyz
            ) != 0)
            {
                goto cleanup;
            }
        }

        internal_circuit_shunt_print_info(
            &shunt
        );

        if (transient_system_create_stage2(
            &stage2,
            &open_model.L,
            &open_model.R,
            &open_model.incidence,
            &open_model.P,
            &shunt,
            transient_scheme_from_config(
                config
            ),
            config->dt
        ) != 0)
        {
            goto cleanup;
        }
    }
    else
    {
        if (transient_system_create_stage1(
            &stage2,
            &open_model.L,
            &open_model.R,
            &open_model.incidence,
            &open_model.P,
            transient_scheme_from_config(
                config
            ),
            config->dt
        ) != 0)
        {
            goto cleanup;
        }
    }

    if (transient_factorize(
        &stage2,
        config->parallel_threads
    ) != 0)
    {
        goto cleanup;
    }

    /*
     * --------------------------------------------------------
     * STAGE-2 A/B PROBE
     * --------------------------------------------------------
     *
     * If a real shunt exists, use its mapped nodes.
     *
     * If --shunt-a-node/--shunt-b-node are supplied WITHOUT
     * --shunt-R, use them only as a voltage probe.  This gives
     * the true open-circuit voltage:
     *
     *     V_oc = phi_A - phi_B.
     */
    if (stage2_has_shunt)
    {
        stage2_probe_enabled =
        1;

        stage2_probe_a =
        shunt.node_a;

        stage2_probe_b =
        shunt.node_b;
    }
    else if (config->shunt_a_node_set &&
        config->shunt_b_node_set)
    {
        if ((size_t)config->shunt_a_node >= open_model.mesh.n_nodes ||
            (size_t)config->shunt_b_node >= open_model.mesh.n_nodes)
        {
            fprintf(
                stderr,
                "ERROR: Stage-2 probe node index is outside open mesh: "
                "A=%d, B=%d, n_nodes=%zu.\n",
                config->shunt_a_node,
                config->shunt_b_node,
                open_model.mesh.n_nodes
            );

            goto cleanup;
        }

        stage2_probe_enabled =
        1;

        stage2_probe_a =
        (size_t)config->shunt_a_node;

        stage2_probe_b =
        (size_t)config->shunt_b_node;

        printf("\n");
        printf("============================================================\n");
        printf("STAGE-2 OPEN-CIRCUIT VOLTAGE PROBE\n");
        printf("============================================================\n");
        printf("A node             : %zu\n", stage2_probe_a);
        printf("B node             : %zu\n", stage2_probe_b);
        printf("Definition         : V_AB = phi_A - phi_B\n");
        printf("Shunt branch       : DISABLED\n");
        printf("Purpose            : measure V_oc(t)\n");
        printf("============================================================\n");
    }


    /*
     * --------------------------------------------------------
     * STATE ARRAYS
     * --------------------------------------------------------
     */
    I1 = calloc(
        closed_model.mesh.n_edges,
        sizeof(double)
    );

    V1 = calloc(
        closed_model.mesh.n_nodes,
        sizeof(double)
    );

    lightning_n = calloc(
        closed_model.mesh.n_nodes,
        sizeof(double)
    );

    lightning_np1 = calloc(
        closed_model.mesh.n_nodes,
        sizeof(double)
    );

    U1_n = calloc(
        closed_model.mesh.n_edges,
        sizeof(double)
    );

    U1_np1 = calloc(
        closed_model.mesh.n_edges,
        sizeof(double)
    );

    I2 = calloc(
        open_model.mesh.n_edges,
        sizeof(double)
    );

    V2 = calloc(
        open_model.mesh.n_nodes,
        sizeof(double)
    );

    I_gamma_n = calloc(
        open_model.mesh.n_nodes,
        sizeof(double)
    );

    I_gamma_np1 = calloc(
        open_model.mesh.n_nodes,
        sizeof(double)
    );

    if (I1 == NULL ||
        V1 == NULL ||
        lightning_n == NULL ||
        lightning_np1 == NULL ||
        U1_n == NULL ||
        U1_np1 == NULL ||
        I2 == NULL ||
        V2 == NULL ||
        I_gamma_n == NULL ||
        I_gamma_np1 == NULL)
    {
        goto cleanup;
    }

    aperture_peak_edge_current = calloc(
        closed_model.mesh.n_edges,
        sizeof(double)
    );

    aperture_peak_gamma = calloc(
        open_model.mesh.n_nodes,
        sizeof(double)
    );

    if (aperture_peak_edge_current == NULL ||
        aperture_peak_gamma == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate aperture diagnostic snapshot.\n"
        );

        goto cleanup;
    }


    /*
     * --------------------------------------------------------
     * TRANSIENT OUTPUT BUFFERS
     * --------------------------------------------------------
     *
     * Внутренняя неизвестная transient:
     *
     *     V_c.
     *
     * Для VTK восстанавливаем:
     *
     *     phi = S V_c,
     *
     *     Q   = F V_c.
     */
    if (config->write_vtk)
    {
        phi1 = calloc(
            closed_model.mesh.n_nodes,
            sizeof(double)
        );

        Q1 = calloc(
            closed_model.mesh.n_nodes,
            sizeof(double)
        );

        Q2 = calloc(
            open_model.mesh.n_nodes,
            sizeof(double)
        );

        if (phi1 == NULL ||
            Q1 == NULL ||
            Q2 == NULL)
        {
            goto cleanup;
        }
    }

    /*
     * Полный потенциал Stage 2 также нужен для проверки
     * внутреннего шунта, даже если VTK выключен.
     */
    if (config->write_vtk ||
        stage2_has_shunt ||
        stage2_probe_enabled)
    {
        phi2 = calloc(
            open_model.mesh.n_nodes,
            sizeof(double)
        );

        if (phi2 == NULL)
        {
            goto cleanup;
        }
    }

    if (config->write_vtk)
    {
        if (write_transient_vtk(
            config,
            "lightning_stage1",
            0,
            0.0,
            &closed_model,
            &stage1,
            I1,
            V1,
            phi1,
            Q1
        ) != 0 ||
        write_transient_vtk(
            config,
            "lightning_stage2",
            0,
            0.0,
            &open_model,
            &stage2,
            I2,
            V2,
            phi2,
            Q2
        ) != 0)
        {
            goto cleanup;
        }

        if (stage2_has_shunt)
        {
            if (write_shunt_visualization_metadata(
                config,
                "lightning_stage2",
                0,
                0.0,
                &open_model.mesh,
                &shunt,
                I_shunt
            ) != 0)
            {
                fprintf(
                    stderr,
                    "ERROR: cannot write Stage-2 shunt visualization metadata.\n"
                );

                goto cleanup;
            }
        }
    }


    /*
     * Initial aperture source:
     *
     *     I1(0) = 0
     *
     * поэтому:
     *
     *     I_gamma(0) = 0.
     */
    if (aperture_coupling_apply(
        &coupling,
        I1,
        I_gamma_n
    ) != 0)
    {
        goto cleanup;
    }


    /*
     * --------------------------------------------------------
     * TIME LOOP
     * --------------------------------------------------------
     */
    size_t n_steps =
    0;

    if (compute_constant_dt_steps(
        config->t_end,
        config->dt,
        &n_steps
    ) != 0)
    {
        goto cleanup;
    }

    size_t frame =
    1;

    size_t progress_step =
    n_steps / 20;

    if (progress_step < 1)
    {
        progress_step =
        1;
    }

    for (size_t step = 0;
         step < n_steps;
    ++step)
         {
             double t_n =
             step
             * config->dt;

             double t_np1 =
             (step + 1)
             * config->dt;

             if (t_np1 >
                 config->t_end)
             {
                 t_np1 =
                 config->t_end;
             }

             if (fabs(
                 t_np1 - t_n - config->dt
             ) >
             1.0e-12
             * fmax(
                 1.0,
                 config->dt
             ))
             {
                 fprintf(
                     stderr,
                     "ERROR: t_end must be an integer multiple of dt.\n"
                 );

                 goto cleanup;
             }


             /*
              * ====================================================
              * STAGE 1
              * ====================================================
              */
             const double *stage1_node_source_n =
             NULL;

             const double *stage1_node_source_np1 =
             NULL;

             const double *stage1_edge_voltage_n =
             NULL;

             const double *stage1_edge_voltage_np1 =
             NULL;

             if (config->excitation ==
                 PEEC_EXCITATION_LIGHTNING_CURRENT)
             {
                 if (lightning_build_node_current(
                     &lightning,
                     closed_model.mesh.n_nodes,
                     t_n,
                     lightning_n
                 ) != 0 ||
                 lightning_build_node_current(
                     &lightning,
                     closed_model.mesh.n_nodes,
                     t_np1,
                     lightning_np1
                 ) != 0)
                 {
                     goto cleanup;
                 }

                 stage1_node_source_n =
                 lightning_n;

                 stage1_node_source_np1 =
                 lightning_np1;
             }
             else
             {
                 if (scattering_build_time_excitation(
                     &closed_model.mesh,
                     &incident_pulse,
                     &scattering_options,
                     t_n,
                     U1_n
                 ) != 0 ||
                 scattering_build_time_excitation(
                     &closed_model.mesh,
                     &incident_pulse,
                     &scattering_options,
                     t_np1,
                     U1_np1
                 ) != 0)
                 {
                     fprintf(
                         stderr,
                         "ERROR: cannot build Stage-1 incident-pulse excitation.\n"
                     );

                     goto cleanup;
                 }

                 stage1_edge_voltage_n =
                 U1_n;

                 stage1_edge_voltage_np1 =
                 U1_np1;
             }

             double strike_Vc_n =
             0.0;

             double return_Vc_n =
             0.0;

             double strike_ATI_n =
             0.0;

             double return_ATI_n =
             0.0;

             if (config->excitation ==
                 PEEC_EXCITATION_LIGHTNING_CURRENT)
             {
                 strike_Vc_n =
                 V1[
                     lightning.strike_node
                 ];

                 return_Vc_n =
                 V1[
                     lightning.return_node
                 ];

                 strike_ATI_n =
                 node_incidence_current(
                     &closed_model.incidence,
                     I1,
                     lightning.strike_node
                 );

                 return_ATI_n =
                 node_incidence_current(
                     &closed_model.incidence,
                     I1,
                     lightning.return_node
                 );
             }

             if (transient_build_rhs_stage1(
                 &stage1,
                 I1,
                 V1,
                 stage1_node_source_n,
                 stage1_node_source_np1,
                 stage1_edge_voltage_n,
                 stage1_edge_voltage_np1
             ) != 0 ||
             transient_solve_step(
                 &stage1
             ) != 0 ||
             transient_extract_stage1_solution(
                 &stage1,
                 I1,
                 V1
             ) != 0)
             {
                 goto cleanup;
             }


             /*
              * ====================================================
              * STAGE-1 KCL DIAGNOSTICS
              * ====================================================
              *
              * Проверяем дискретное узловое уравнение:
              *
              *     A^T I - F dV_c/dt + J = 0.
              *
              * Для Backward Euler:
              *
              *     I_theta = I^(n+1),
              *     J_theta = J^(n+1).
              *
              * Для trapezoidal:
              *
              *     I_theta = (I^n + I^(n+1))/2,
              *     J_theta = (J^n + J^(n+1))/2.
              */
             if (config->excitation ==
                 PEEC_EXCITATION_LIGHTNING_CURRENT)
             {
                 double theta =
                 config->time_order == 1
                 ? 1.0
                 : 0.5;

                 size_t strike =
                 lightning.strike_node;

                 size_t return_node =
                 lightning.return_node;

                 double strike_ATI_np1 =
                 node_incidence_current(
                     &closed_model.incidence,
                     I1,
                     strike
                 );

                 double return_ATI_np1 =
                 node_incidence_current(
                     &closed_model.incidence,
                     I1,
                     return_node
                 );

                 double strike_ATI_theta =
                 (1.0 - theta)
                 * strike_ATI_n
                 +
                 theta
                 * strike_ATI_np1;

                 double return_ATI_theta =
                 (1.0 - theta)
                 * return_ATI_n
                 +
                 theta
                 * return_ATI_np1;

                 double strike_J_theta =
                 (1.0 - theta)
                 * lightning_n[strike]
                 +
                 theta
                 * lightning_np1[strike];

                 double return_J_theta =
                 (1.0 - theta)
                 * lightning_n[return_node]
                 +
                 theta
                 * lightning_np1[return_node];

                 double F_strike =
                 1.0
                 /
                 closed_model.P.data[
                     strike
                     * closed_model.P.n
                     + strike
                 ];

                 double F_return =
                 1.0
                 /
                 closed_model.P.data[
                     return_node
                     * closed_model.P.n
                     + return_node
                 ];

                 double strike_Cdot =
                 F_strike
                 * (
                     V1[strike]
                     - strike_Vc_n
                 )
                 / config->dt;

                 double return_Cdot =
                 F_return
                 * (
                     V1[return_node]
                     - return_Vc_n
                 )
                 / config->dt;

                 double strike_residual =
                 strike_ATI_theta
                 - strike_Cdot
                 + strike_J_theta;

                 double return_residual =
                 return_ATI_theta
                 - return_Cdot
                 + return_J_theta;

                 double source_abs =
                 fabs(
                     lightning_current(
                         &lightning,
                         t_np1
                     )
                 );

                 if (source_abs >
                     sampled_source_max)
                 {
                     sampled_source_max =
                     source_abs;

                     sampled_source_time =
                     t_np1;

                     sampled_strike_ATI =
                     strike_ATI_theta;

                     sampled_strike_Cdot =
                     strike_Cdot;

                     sampled_strike_J =
                     strike_J_theta;

                     sampled_strike_residual =
                     strike_residual;

                     sampled_return_ATI =
                     return_ATI_theta;

                     sampled_return_Cdot =
                     return_Cdot;

                     sampled_return_J =
                     return_J_theta;

                     sampled_return_residual =
                     return_residual;
                 }

                 size_t max_edge =
                 0;

                 double max_edge_current =
                 edge_current_max_abs(
                     I1,
                     closed_model.mesh.n_edges,
                     &max_edge
                 );

                 if (max_edge_current >
                     stage1_edge_max)
                 {
                     stage1_edge_max =
                     max_edge_current;

                     stage1_edge_max_time =
                     t_np1;

                     stage1_edge_max_index =
                     max_edge;
                 }
             }


             /*
              * I_gamma^(n+1)
              */
             if (aperture_coupling_apply(
                 &coupling,
                 I1,
                 I_gamma_np1
             ) != 0)
             {
                 goto cleanup;
             }

             /*
              * sampled_source_time is updated above whenever the
              * discrete lightning-source maximum increases.
              * Therefore this snapshot is overwritten until the
              * final sampled source peak is reached.
              */
             if (config->excitation ==
                 PEEC_EXCITATION_LIGHTNING_CURRENT &&
                 fabs(
                     t_np1
                     -
                     sampled_source_time
                 )
                 <=
                 0.25
                 * config->dt)
             {
                 memcpy(
                     aperture_peak_edge_current,
                     I1,
                     closed_model.mesh.n_edges
                     * sizeof(double)
                 );

                 memcpy(
                     aperture_peak_gamma,
                     I_gamma_np1,
                     open_model.mesh.n_nodes
                     * sizeof(double)
                 );

                 aperture_peak_time =
                 t_np1;
             }


             /*
              * ====================================================
              * STAGE 2
              * ====================================================
              */
             if (stage2_has_shunt)
             {
                 if (transient_build_rhs_stage2(
                     &stage2,
                     I2,
                     V2,
                     I_shunt,
                     I_gamma_n,
                     I_gamma_np1,
                     NULL,
                     NULL
                 ) != 0 ||
                 transient_solve_step(
                     &stage2
                 ) != 0 ||
                 transient_extract_stage2_solution(
                     &stage2,
                     I2,
                     V2,
                     &I_shunt
                 ) != 0)
                 {
                     goto cleanup;
                 }
             }
             else
             {
                 if (transient_build_rhs_stage1(
                     &stage2,
                     I2,
                     V2,
                     I_gamma_n,
                     I_gamma_np1,
                     NULL,
                     NULL
                 ) != 0 ||
                 transient_solve_step(
                     &stage2
                 ) != 0 ||
                 transient_extract_stage1_solution(
                     &stage2,
                     I2,
                     V2
                 ) != 0)
                 {
                     goto cleanup;
                 }
             }


             /*
              * ====================================================
              * STAGE-2 APERTURE / PORT DIAGNOSTICS
              * ====================================================
              */
             {
                 double gamma_max_step =
                 0.0;

                 size_t gamma_max_node_step =
                 0;

                 double gamma_abs_sum =
                 0.0;

                 for (size_t j = 0;
                      j < open_model.mesh.n_nodes;
                 ++j)
                      {
                          double value =
                          fabs(
                              I_gamma_np1[j]
                          );

                          gamma_abs_sum +=
                          value;

                          if (value >
                              gamma_max_step)
                          {
                              gamma_max_step =
                              value;

                              gamma_max_node_step =
                              j;
                          }
                      }

                      /*
                       * For a balanced injection/extraction pattern,
                       *
                       *     0.5 * sum |I_gamma,j|
                       *
                       * is a useful scalar measure of transferred
                       * current magnitude; unlike sum(I_gamma), it
                       * does not cancel to approximately zero.
                       */
                      double gamma_half_l1 =
                      0.5
                      * gamma_abs_sum;

                      if (gamma_max_step >
                          gamma_max_global)
                      {
                          gamma_max_global =
                          gamma_max_step;

                          gamma_max_global_time =
                          t_np1;

                          gamma_max_global_node =
                          gamma_max_node_step;

                          gamma_half_l1_at_global_max =
                          gamma_half_l1;
                      }

                      if (gamma_half_l1 >
                          gamma_half_l1_global)
                      {
                          gamma_half_l1_global =
                          gamma_half_l1;

                          gamma_half_l1_global_time =
                          t_np1;
                      }

                      if (stage2_probe_enabled)
                      {
                          if (transient_recover_node_state(
                              &stage2,
                              V2,
                              phi2,
                              NULL
                          ) != 0)
                          {
                              goto cleanup;
                          }

                          double port_voltage =
                          phi2[stage2_probe_a]
                          -
                          phi2[stage2_probe_b];

                          if (fabs(port_voltage) >
                              port_voltage_abs_max)
                          {
                              port_voltage_abs_max =
                              fabs(
                                  port_voltage
                              );

                              port_voltage_at_abs_max =
                              port_voltage;

                              port_voltage_abs_max_time =
                              t_np1;
                          }
                      }

                      if (stage2_has_shunt &&
                          fabs(I_shunt) >
                          shunt_current_abs_max)
                      {
                          shunt_current_abs_max =
                          fabs(
                              I_shunt
                          );

                          shunt_current_at_abs_max =
                          I_shunt;

                          shunt_current_abs_max_time =
                          t_np1;
                      }
             }


             /*
              * Следующий шаг:
              *
              *     I_gamma^n <- I_gamma^(n+1).
              */
             memcpy(
                 I_gamma_n,
                 I_gamma_np1,
                 open_model.mesh.n_nodes
                 * sizeof(double)
             );


             /*
              * ====================================================
              * OUTPUT
              * ====================================================
              */
             if (config->write_vtk &&
                 (
                     (step + 1) % (size_t)config->vtk_every == 0 ||
                     step + 1 == n_steps
                 ))
             {
                 if (write_transient_vtk(
                     config,
                     "lightning_stage1",
                     frame,
                     t_np1,
                     &closed_model,
                     &stage1,
                     I1,
                     V1,
                     phi1,
                     Q1
                 ) != 0 ||
                 write_transient_vtk(
                     config,
                     "lightning_stage2",
                     frame,
                     t_np1,
                     &open_model,
                     &stage2,
                     I2,
                     V2,
                     phi2,
                     Q2
                 ) != 0)
                 {
                     goto cleanup;
                 }

                 if (stage2_has_shunt)
                 {
                     if (write_shunt_visualization_metadata(
                         config,
                         "lightning_stage2",
                         frame,
                         t_np1,
                         &open_model.mesh,
                         &shunt,
                         I_shunt
                     ) != 0)
                     {
                         fprintf(
                             stderr,
                             "ERROR: cannot write Stage-2 shunt visualization metadata.\n"
                         );

                         goto cleanup;
                     }
                 }

                 ++frame;
             }


             /*
              * ====================================================
              * PROGRESS
              * ====================================================
              */
             if ((step + 1) % progress_step == 0 ||
                 step + 1 == n_steps)
             {
                 double net_gamma =
                 0.0;

                 double max_abs_gamma =
                 0.0;

                 double half_l1_gamma =
                 0.0;

                 for (size_t j = 0;
                      j < open_model.mesh.n_nodes;
                 ++j)
                      {
                          double gamma =
                          I_gamma_np1[j];

                          double gamma_abs =
                          fabs(
                              gamma
                          );

                          net_gamma +=
                          gamma;

                          half_l1_gamma +=
                          gamma_abs;

                          if (gamma_abs >
                              max_abs_gamma)
                          {
                              max_abs_gamma =
                              gamma_abs;
                          }
                      }

                      half_l1_gamma *=
                      0.5;

                      if (config->excitation ==
                          PEEC_EXCITATION_LIGHTNING_CURRENT)
                      {
                          printf(
                              "Transient 2-stage: %6.1f%%  "
                              "t=%.9e s  "
                              "I_L=%.9e A  "
                              "net I_gamma=%.9e A  "
                              "max|I_gamma|=%.9e A  "
                              "0.5sum|I_gamma|=%.9e A",
                              100.0
                              * (double)(step + 1)
                              / (double)n_steps,
                                 t_np1,
                                 lightning_current(
                                     &lightning,
                                     t_np1
                                 ),
                                 net_gamma,
                                 max_abs_gamma,
                                 half_l1_gamma
                          );
                      }
                      else
                      {
                          printf(
                              "Transient 2-stage: %6.1f%%  "
                              "t=%.9e s  "
                              "incident-pulse  "
                              "net I_gamma=%.9e A  "
                              "max|I_gamma|=%.9e A  "
                              "0.5sum|I_gamma|=%.9e A",
                              100.0
                              * (double)(step + 1)
                              / (double)n_steps,
                                 t_np1,
                                 net_gamma,
                                 max_abs_gamma,
                                 half_l1_gamma
                          );
                      }

                      if (stage2_probe_enabled)
                      {
                          if (transient_recover_node_state(
                              &stage2,
                              V2,
                              phi2,
                              NULL
                          ) != 0)
                          {
                              goto cleanup;
                          }

                          double port_voltage =
                          phi2[stage2_probe_a]
                          -
                          phi2[stage2_probe_b];

                          printf(
                              "  V_AB=%.9e V",
                              port_voltage
                          );
                      }

                      if (stage2_has_shunt)
                      {
                          double shunt_voltage =
                          internal_circuit_shunt_voltage(
                              &shunt,
                              phi2,
                              open_model.mesh.n_nodes
                          );

                          double shunt_error =
                          internal_circuit_shunt_residual(
                              &shunt,
                              phi2,
                              open_model.mesh.n_nodes,
                              I_shunt
                          );

                          printf(
                              "  I_shunt=%.9e A  "
                              "V_shunt=%.9e V  "
                              "Ohm residual=%.3e V",
                              I_shunt,
                              shunt_voltage,
                              shunt_error
                          );
                      }

                      printf("\n");
             }
         }

         if (aperture_peak_edge_current != NULL &&
             aperture_peak_gamma != NULL)
         {
             print_aperture_current_transfer_diagnostic(
                 &coupling,
                 &closed_model.mesh,
                 aperture_peak_edge_current,
                 aperture_peak_gamma,
                 aperture_peak_time
             );
         }

         printf("\n");
         printf("============================================================\n");
         printf("STAGE-2 APERTURE / PORT DIAGNOSTICS\n");
         printf("============================================================\n");

         printf(
             "Global max |I_gamma,j|    : %.9e A\n",
             gamma_max_global
         );

         printf(
             "I_gamma max node          : %zu\n",
             gamma_max_global_node
         );

         printf(
             "I_gamma max time          : %.9e s (%.6f us)\n",
                gamma_max_global_time,
                gamma_max_global_time * 1.0e6
         );

         printf(
             "0.5 sum|I_gamma| there    : %.9e A\n",
             gamma_half_l1_at_global_max
         );

         printf(
             "Global max 0.5sum|I_gamma|: %.9e A\n",
             gamma_half_l1_global
         );

         printf(
             "0.5sum max time           : %.9e s (%.6f us)\n",
                gamma_half_l1_global_time,
                gamma_half_l1_global_time * 1.0e6
         );

         if (stage2_probe_enabled)
         {
             printf("\n");
             printf(
                 "Port A/B                 : %zu / %zu\n",
                 stage2_probe_a,
                 stage2_probe_b
             );

             printf(
                 "max |V_AB|               : %.9e V\n",
                 port_voltage_abs_max
             );

             printf(
                 "V_AB at that time        : %.9e V\n",
                 port_voltage_at_abs_max
             );

             printf(
                 "V_AB max time            : %.9e s (%.6f us)\n",
                    port_voltage_abs_max_time,
                    port_voltage_abs_max_time * 1.0e6
             );

             if (!stage2_has_shunt)
             {
                 printf(
                     "Port state               : OPEN CIRCUIT (V_AB = V_oc)\n"
                 );
             }
         }

         if (stage2_has_shunt)
         {
             printf("\n");

             printf(
                 "Shunt resistance         : %.9e Ohm\n",
                 shunt.resistance
             );

             printf(
                 "max |I_shunt|            : %.9e A\n",
                 shunt_current_abs_max
             );

             printf(
                 "I_shunt at that time     : %.9e A\n",
                 shunt_current_at_abs_max
             );

             printf(
                 "I_shunt max time         : %.9e s (%.6f us)\n",
                    shunt_current_abs_max_time,
                    shunt_current_abs_max_time * 1.0e6
             );

             if (shunt.resistance > 0.0)
             {
                 printf(
                     "R * I_shunt(max-I time)  : %.9e V\n",
                        shunt.resistance
                        * shunt_current_at_abs_max
                 );
             }
         }

         printf("============================================================\n");

         if (config->excitation ==
             PEEC_EXCITATION_LIGHTNING_CURRENT)
         {
             printf("\n");
             printf("============================================================\n");
             printf("STAGE-1 LIGHTNING / KCL DIAGNOSTICS\n");
             printf("============================================================\n");

             printf(
                 "Peak sampled source time : %.9e s (%.6f us)\n",
                    sampled_source_time,
                    sampled_source_time * 1.0e6
             );

             printf(
                 "Peak sampled |I_L|       : %.9e A\n",
                 sampled_source_max
             );

             printf("\nStrike node %zu\n",
                    lightning.strike_node);

             printf(
                 "  J_source(theta)        : %.9e A\n",
                    sampled_strike_J
             );

             printf(
                 "  (A^T I)(theta)         : %.9e A\n",
                    sampled_strike_ATI
             );

             printf(
                 "  F dV_c/dt              : %.9e A\n",
                 sampled_strike_Cdot
             );

             printf(
                 "  KCL residual           : %.9e A\n",
                 sampled_strike_residual
             );

             printf("\nReturn node %zu\n",
                    lightning.return_node);

             printf(
                 "  J_source(theta)        : %.9e A\n",
                    sampled_return_J
             );

             printf(
                 "  (A^T I)(theta)         : %.9e A\n",
                    sampled_return_ATI
             );

             printf(
                 "  F dV_c/dt              : %.9e A\n",
                 sampled_return_Cdot
             );

             printf(
                 "  KCL residual           : %.9e A\n",
                 sampled_return_residual
             );

             printf("\n");

             printf(
                 "Global max |edge current|: %.9e A\n",
                 stage1_edge_max
             );

             printf(
                 "Max edge index           : %zu\n",
                 stage1_edge_max_index
             );

             printf(
                 "Max edge-current time    : %.9e s (%.6f us)\n",
                    stage1_edge_max_time,
                    stage1_edge_max_time * 1.0e6
             );

             printf("============================================================\n");
         }

         result =
         0;


         cleanup:

         free(
             U1_np1
         );

         free(
             U1_n
         );

         free(
             aperture_peak_edge_current
         );

         free(
             aperture_peak_gamma
         );

         free(
             I_gamma_np1
         );

         free(
             I_gamma_n
         );

         free(
             Q2
         );

         free(
             phi2
         );

         free(
             V2
         );

         free(
             I2
         );

         free(
             lightning_np1
         );

         free(
             lightning_n
         );

         free(
             Q1
         );

         free(
             phi1
         );

         free(
             V1
         );

         free(
             I1
         );

         transient_system_free(
             &stage2
         );

         transient_system_free(
             &stage1
         );

         aperture_coupling_free(
             &coupling
         );

         peec_model_free(
             &open_model
         );

         peec_model_free(
             &closed_model
         );

         return result;
}


/*
 * ============================================================
 * LIGHTNING / ONE-STAGE DISTRIBUTED SLOT CELLS
 * ============================================================
 *
 * One open PEEC mesh is solved.
 *
 * The aperture is replaced by N local parallel LC cells:
 *
 *     (A_m) ---- L_m || C_m ---- (B_m)
 *
 * The state is:
 *
 *     [ I_edge, V_c, I_L,slot, optional I_shunt ].
 *
 * No Stage-1 -> Stage-2 current transfer is used.
 */
static int run_lightning_slot_cells(
    const PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    PeecModel model;

    peec_model_init(
        &model
    );

    SlotCircuit slot;

    slot_circuit_init(
        &slot
    );

    SlotLineModel slot_line;

    memset(
        &slot_line,
        0,
        sizeof(slot_line)
    );

    InternalCircuitShunt shunt;

    internal_circuit_shunt_default(
        &shunt
    );

    TransientSystem transient;

    transient_system_init(
        &transient
    );

    LightningSource lightning;

    IncidentField incident_pulse;

    incident_field_init(
        &incident_pulse
    );

    ScatteringOptions scattering_options;

    scattering_options_default(
        &scattering_options,
        config->parallel_threads
    );

    double *I =
    NULL;

    double *V =
    NULL;

    double *I_slot =
    NULL;

    double *phi =
    NULL;

    double *Q =
    NULL;

    double *source_n =
    NULL;

    double *source_np1 =
    NULL;

    double *edge_voltage_n =
    NULL;

    double *edge_voltage_np1 =
    NULL;

    double I_shunt =
    0.0;

    int has_shunt =
    config->use_internal_shunt;

    int result =
    -1;

    double max_slot_current =
    0.0;

    size_t max_slot_cell =
    0;

    double max_slot_time =
    0.0;

    double max_shunt_current =
    0.0;

    double max_shunt_time =
    0.0;

    printf("\n");
    printf("============================================================\n");
    printf("LIGHTNING / ONE-STAGE DISTRIBUTED SLOT-CELL MODEL\n");
    printf("============================================================\n");

    /*
     * Only the physical open-body mesh is assembled.
     */
    if (peec_assemble_model(
        config->open_mesh_file,
        config->parallel_threads,
        &model
    ) != 0)
    {
        goto cleanup;
    }

    if (peec_export_matrices_if_requested(
        config,
        &model,
        config->matrix_directory
    ) != 0)
    {
        goto cleanup;
    }

    /*
     * Equivalent aperture line.
     */
    if (slot_line_model_build(
        &slot_line,
        config->slot_width,
        config->slot_wall_span,
        config->slot_wall_thickness,
        EPSILON0 * config->slot_epsilon_r,
        MU0 * config->slot_mu_r
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot build slot transmission-line model.\n"
        );

        goto cleanup;
    }

    if (slot_circuit_load_cells(
        &slot,
        config->slot_cells_file,
        &slot_line,
        model.mesh.n_nodes
    ) != 0)
    {
        goto cleanup;
    }

    slot_circuit_print_info(
        &slot
    );

    /*
     * External excitation is applied directly to the one-stage
     * open PEEC model.
     */
    if (config->excitation ==
        PEEC_EXCITATION_LIGHTNING_CURRENT)
    {
        if (configure_lightning_source(
            config,
            &model.mesh,
            &lightning
        ) != 0)
        {
            goto cleanup;
        }
    }
    else
    {
        if (configure_incident_pulse(
            config,
            &incident_pulse
        ) != 0)
        {
            goto cleanup;
        }
    }

    /*
     * Optional internal resistive load.
     */
    const InternalCircuitShunt *shunt_ptr =
    NULL;

    if (has_shunt)
    {
        shunt.resistance =
        config->shunt_resistance;

        if (config->shunt_a_node_set &&
            config->shunt_b_node_set)
        {
            if ((size_t)config->shunt_a_node >= model.mesh.n_nodes ||
                (size_t)config->shunt_b_node >= model.mesh.n_nodes)
            {
                fprintf(
                    stderr,
                    "ERROR: shunt node index is outside slot-cell mesh.\n"
                );

                goto cleanup;
            }

            if (internal_circuit_shunt_set_nodes(
                &shunt,
                &model.mesh,
                (size_t)config->shunt_a_node,
                                                 (size_t)config->shunt_b_node
            ) != 0)
            {
                goto cleanup;
            }
        }
        else
        {
            if (internal_circuit_shunt_map_coordinates(
                &shunt,
                &model.mesh,
                config->shunt_a_xyz,
                config->shunt_b_xyz
            ) != 0)
            {
                goto cleanup;
            }
        }

        internal_circuit_shunt_print_info(
            &shunt
        );

        shunt_ptr =
        &shunt;
    }

    if (transient_system_create_slot_cells(
        &transient,
        &model.L,
        &model.R,
        &model.incidence,
        &model.P,
        &slot,
        shunt_ptr,
        transient_scheme_from_config(
            config
        ),
        config->dt
    ) != 0)
    {
        goto cleanup;
    }

    if (transient_factorize(
        &transient,
        config->parallel_threads
    ) != 0)
    {
        goto cleanup;
    }

    transient_print_info(
        &transient
    );

    I = calloc(
        model.mesh.n_edges,
        sizeof(double)
    );

    V = calloc(
        model.mesh.n_nodes,
        sizeof(double)
    );

    I_slot = calloc(
        slot.n_cells,
        sizeof(double)
    );

    source_n = calloc(
        model.mesh.n_nodes,
        sizeof(double)
    );

    source_np1 = calloc(
        model.mesh.n_nodes,
        sizeof(double)
    );

    edge_voltage_n = calloc(
        model.mesh.n_edges,
        sizeof(double)
    );

    edge_voltage_np1 = calloc(
        model.mesh.n_edges,
        sizeof(double)
    );

    if (I == NULL ||
        V == NULL ||
        I_slot == NULL ||
        source_n == NULL ||
        source_np1 == NULL ||
        edge_voltage_n == NULL ||
        edge_voltage_np1 == NULL)
    {
        goto cleanup;
    }

    if (config->write_vtk ||
        has_shunt)
    {
        phi = calloc(
            model.mesh.n_nodes,
            sizeof(double)
        );

        if (phi == NULL)
        {
            goto cleanup;
        }
    }

    if (config->write_vtk)
    {
        Q = calloc(
            model.mesh.n_nodes,
            sizeof(double)
        );

        if (Q == NULL)
        {
            goto cleanup;
        }

        if (write_transient_vtk(
            config,
            "lightning_slot_cells",
            0,
            0.0,
            &model,
            &transient,
            I,
            V,
            phi,
            Q
        ) != 0 ||
        write_slot_visualization_metadata(
            config,
            "lightning_slot_cells",
            0,
            0.0,
            &slot,
            I_slot
        ) != 0)
        {
            goto cleanup;
        }

        if (has_shunt)
        {
            if (write_shunt_visualization_metadata(
                config,
                "lightning_slot_cells",
                0,
                0.0,
                &model.mesh,
                &shunt,
                I_shunt
            ) != 0)
            {
                goto cleanup;
            }
        }
    }

    size_t n_steps =
    0;

    if (compute_constant_dt_steps(
        config->t_end,
        config->dt,
        &n_steps
    ) != 0)
    {
        goto cleanup;
    }

    size_t frame =
    1;

    size_t progress_step =
    n_steps / 20;

    if (progress_step < 1)
    {
        progress_step =
        1;
    }

    for (size_t step = 0;
         step < n_steps;
    ++step)
         {
             double t_n =
             step
             * config->dt;

             double t_np1 =
             (step + 1)
             * config->dt;

             const double *node_source_n_ptr =
             NULL;

             const double *node_source_np1_ptr =
             NULL;

             const double *edge_voltage_n_ptr =
             NULL;

             const double *edge_voltage_np1_ptr =
             NULL;

             if (config->excitation ==
                 PEEC_EXCITATION_LIGHTNING_CURRENT)
             {
                 if (lightning_build_node_current(
                     &lightning,
                     model.mesh.n_nodes,
                     t_n,
                     source_n
                 ) != 0 ||
                 lightning_build_node_current(
                     &lightning,
                     model.mesh.n_nodes,
                     t_np1,
                     source_np1
                 ) != 0)
                 {
                     goto cleanup;
                 }

                 node_source_n_ptr =
                 source_n;

                 node_source_np1_ptr =
                 source_np1;
             }
             else
             {
                 if (scattering_build_time_excitation(
                     &model.mesh,
                     &incident_pulse,
                     &scattering_options,
                     t_n,
                     edge_voltage_n
                 ) != 0 ||
                 scattering_build_time_excitation(
                     &model.mesh,
                     &incident_pulse,
                     &scattering_options,
                     t_np1,
                     edge_voltage_np1
                 ) != 0)
                 {
                     goto cleanup;
                 }

                 edge_voltage_n_ptr =
                 edge_voltage_n;

                 edge_voltage_np1_ptr =
                 edge_voltage_np1;
             }

             if (transient_build_rhs_slot_cells(
                 &transient,
                 I,
                 V,
                 I_slot,
                 I_shunt,
                 node_source_n_ptr,
                 node_source_np1_ptr,
                 edge_voltage_n_ptr,
                 edge_voltage_np1_ptr
             ) != 0 ||
             transient_solve_step(
                 &transient
             ) != 0 ||
             transient_extract_slot_cells_solution(
                 &transient,
                 I,
                 V,
                 I_slot,
                 has_shunt
                 ? &I_shunt
                 : NULL
             ) != 0)
             {
                 goto cleanup;
             }

             for (size_t m = 0;
                  m < slot.n_cells;
             ++m)
                  {
                      if (fabs(I_slot[m]) >
                          max_slot_current)
                      {
                          max_slot_current =
                          fabs(
                              I_slot[m]
                          );

                          max_slot_cell =
                          m;

                          max_slot_time =
                          t_np1;
                      }
                  }

                  if (has_shunt &&
                      fabs(I_shunt) >
                      max_shunt_current)
                  {
                      max_shunt_current =
                      fabs(
                          I_shunt
                      );

                      max_shunt_time =
                      t_np1;
                  }

                  if (config->write_vtk &&
                      (
                          (step + 1)
                          % (size_t)config->vtk_every
                          == 0
                          ||
                          step + 1 == n_steps
                      ))
                  {
                      if (write_transient_vtk(
                          config,
                          "lightning_slot_cells",
                          frame,
                          t_np1,
                          &model,
                          &transient,
                          I,
                          V,
                          phi,
                          Q
                      ) != 0 ||
                      write_slot_visualization_metadata(
                          config,
                          "lightning_slot_cells",
                          frame,
                          t_np1,
                          &slot,
                          I_slot
                      ) != 0)
                      {
                          goto cleanup;
                      }

                      if (has_shunt)
                      {
                          if (write_shunt_visualization_metadata(
                              config,
                              "lightning_slot_cells",
                              frame,
                              t_np1,
                              &model.mesh,
                              &shunt,
                              I_shunt
                          ) != 0)
                          {
                              goto cleanup;
                          }
                      }

                      ++frame;
                  }

                  if ((step + 1) % progress_step == 0 ||
                      step + 1 == n_steps)
                  {
                      double step_slot_max =
                      0.0;

                      for (size_t m = 0;
                           m < slot.n_cells;
                      ++m)
                           {
                               if (fabs(I_slot[m]) >
                                   step_slot_max)
                               {
                                   step_slot_max =
                                   fabs(
                                       I_slot[m]
                                   );
                               }
                           }

                           if (config->excitation ==
                               PEEC_EXCITATION_LIGHTNING_CURRENT)
                           {
                               printf(
                                   "Transient slot-cells: %6.1f%%  "
                                   "t=%.9e s  "
                                   "I_L=%.9e A  "
                                   "max|I_slot|=%.9e A",
                                   100.0
                                   * (double)(step + 1)
                                   / (double)n_steps,
                                      t_np1,
                                      lightning_current(
                                          &lightning,
                                          t_np1
                                      ),
                                      step_slot_max
                               );
                           }
                           else
                           {
                               printf(
                                   "Transient slot-cells: %6.1f%%  "
                                   "t=%.9e s  "
                                   "max|I_slot|=%.9e A",
                                   100.0
                                   * (double)(step + 1)
                                   / (double)n_steps,
                                      t_np1,
                                      step_slot_max
                               );
                           }

                           if (has_shunt)
                           {
                               if (transient_recover_node_state(
                                   &transient,
                                   V,
                                   phi,
                                   NULL
                               ) != 0)
                               {
                                   goto cleanup;
                               }

                               double shunt_voltage =
                               internal_circuit_shunt_voltage(
                                   &shunt,
                                   phi,
                                   model.mesh.n_nodes
                               );

                               printf(
                                   "  I_shunt=%.9e A  V_shunt=%.9e V",
                                   I_shunt,
                                   shunt_voltage
                               );
                           }

                           printf("\n");
                  }
         }

         printf("\n");
         printf("============================================================\n");
         printf("ONE-STAGE SLOT-CELL DIAGNOSTICS\n");
         printf("============================================================\n");

         printf(
             "Maximum |I_slot|   : %.9e A\n",
             max_slot_current
         );

         printf(
             "Slot cell          : %zu\n",
             max_slot_cell
         );

         printf(
             "Time               : %.9e s (%.6f us)\n",
                max_slot_time,
                max_slot_time * 1.0e6
         );

         if (has_shunt)
         {
             printf(
                 "Maximum |I_shunt|  : %.9e A\n",
                 max_shunt_current
             );

             printf(
                 "I_shunt max time   : %.9e s (%.6f us)\n",
                    max_shunt_time,
                    max_shunt_time * 1.0e6
             );
         }

         printf("============================================================\n");

         result =
         0;


         cleanup:

         free(
             edge_voltage_np1
         );

         free(
             edge_voltage_n
         );

         free(
             source_np1
         );

         free(
             source_n
         );

         free(
             Q
         );

         free(
             phi
         );

         free(
             I_slot
         );

         free(
             V
         );

         free(
             I
         );

         transient_system_free(
             &transient
         );

         slot_circuit_free(
             &slot
         );

         peec_model_free(
             &model
         );

         return result;
}


/*
 * ============================================================
 * LIGHTNING DISPATCH
 * ============================================================
 */
static int run_lightning(
    const PeecConfig *config
)
{
    if (config == NULL)
    {
        return -1;
    }

    if (config->physics !=
        PEEC_PHYSICS_QUASISTATIC)
    {
        fprintf(
            stderr,
            "ERROR: transient lightning is currently implemented "
            "for quasistatic PEEC only.\n"
        );

        return -1;
    }

    if (config->lightning_case ==
        PEEC_LIGHTNING_CLOSED_BODY)
    {
        return run_lightning_closed(
            config
        );
    }

    if (config->lightning_case ==
        PEEC_LIGHTNING_TWO_STAGE)
    {
        return run_lightning_two_stage(
            config
        );
    }

    if (config->lightning_case ==
        PEEC_LIGHTNING_SLOT_CELLS)
    {
        return run_lightning_slot_cells(
            config
        );
    }

    return -1;
}



/*
 * ============================================================
 * CONFIG INPUT MODE
 * ============================================================
 *
 * Choose here how the program receives run parameters.
 *
 * PEEC_INPUT_COMMAND_LINE:
 *
 *     ./PeecSolver --task lightning ...
 *
 * PEEC_INPUT_FILE:
 *
 *     PeecSolver.exe
 *
 * and all options are read from PEEC_CONFIG_FILE.
 *
 * The configuration file uses exactly the same option names
 * as the command line, for example:
 *
 *     --task lightning
 *     --dt 1e-7
 *     --vtk
 *
 * Empty lines and comments beginning with # are ignored.
 * Single and double quoted values are supported, which is useful
 * for Windows paths containing spaces.
 */
typedef enum
{
    PEEC_INPUT_COMMAND_LINE = 0,
    PEEC_INPUT_FILE = 1

} PeecInputMode;


/*
 * ============================================================
 * USER-SELECTION AREA
 * ============================================================
 *
 * Change only these two values when switching between CLI and
 * configuration-file operation.
 */
/*
 * Automatic platform-dependent input mode:
 *
 * Windows:
 *     read parameters from PEEC_CONFIG_FILE.
 *
 * Linux / macOS / Unix:
 *     use command-line arguments.
 *
 * The OS is detected by the compiler at build time.
 */
#ifdef _WIN32

static const PeecInputMode PEEC_INPUT_MODE =
PEEC_INPUT_FILE;

#else

static const PeecInputMode PEEC_INPUT_MODE =
PEEC_INPUT_COMMAND_LINE;

#endif


static const char *PEEC_CONFIG_FILE =
"cases/run_lightning.cfg";


         /*
          * ============================================================
          * CONFIG-FILE TOKEN STORAGE
          * ============================================================
          */
         typedef struct
         {
             int argc;
             int capacity;
             char **argv;

         } ConfigArguments;


         static void config_arguments_init(
             ConfigArguments *arguments
         )
         {
             if (arguments == NULL)
             {
                 return;
             }

             arguments->argc =
             0;

             arguments->capacity =
             0;

             arguments->argv =
             NULL;
         }


         static void config_arguments_free(
             ConfigArguments *arguments
         )
         {
             if (arguments == NULL)
             {
                 return;
             }

             for (int i = 0;
                  i < arguments->argc;
             ++i)
                  {
                      free(
                          arguments->argv[i]
                      );
                  }

                  free(
                      arguments->argv
                  );

                  config_arguments_init(
                      arguments
                  );
         }


         static int config_arguments_push(
             ConfigArguments *arguments,
             const char *text,
             size_t length
         )
         {
             if (arguments == NULL ||
                 text == NULL)
             {
                 return -1;
             }

             if (arguments->argc >=
                 arguments->capacity)
             {
                 int new_capacity =
                 arguments->capacity == 0
                 ? 32
                 : 2 * arguments->capacity;

                 char **new_argv =
                 realloc(
                     arguments->argv,
                     (size_t)new_capacity
                     * sizeof(char *)
                 );

                 if (new_argv == NULL)
                 {
                     return -1;
                 }

                 arguments->argv =
                 new_argv;

                 arguments->capacity =
                 new_capacity;
             }

             char *copy =
             malloc(
                 length + 1
             );

             if (copy == NULL)
             {
                 return -1;
             }

             memcpy(
                 copy,
                 text,
                 length
             );

             copy[length] =
             '\0';

             arguments->argv[
                 arguments->argc
             ] =
             copy;

             ++arguments->argc;

             return 0;
         }


         /*
          * ============================================================
          * CONFIG-FILE PARSER
          * ============================================================
          *
          * The file is converted into an artificial argc/argv array and
          * then passed to the existing config_parse_cli().
          *
          * Therefore the interpretation and validation of PEEC options
          * stay in one place: config.c.
          *
          * Supported syntax:
          *
          *     --task lightning
          *     --dt 1e-7
          *     --vtk
          *
          *     # comment
          *
          *     --closed-mesh "C:/PEEC project/meshes/body_closed.msh"
          *
          * Windows backslashes are not treated as escape characters.
          */
         static int config_parse_file_as_cli(
             PeecConfig *config,
             const char *program_name,
             const char *filename
         )
         {
             if (config == NULL ||
                 program_name == NULL ||
                 filename == NULL)
             {
                 return -1;
             }

             FILE *file =
             fopen(
                 filename,
                 "rb"
             );

             if (file == NULL)
             {
                 fprintf(
                     stderr,
                     "ERROR: cannot open configuration file '%s': %s\n",
                     filename,
                     strerror(errno)
                 );

                 return -1;
             }

             if (fseek(
                 file,
                 0,
                 SEEK_END
             ) != 0)
             {
                 fclose(
                     file
                 );

                 return -1;
             }

             long file_size_long =
             ftell(
                 file
             );

             if (file_size_long < 0)
             {
                 fclose(
                     file
                 );

                 return -1;
             }

             if (fseek(
                 file,
                 0,
                 SEEK_SET
             ) != 0)
             {
                 fclose(
                     file
                 );

                 return -1;
             }

             size_t file_size =
             (size_t)file_size_long;

             char *buffer =
             malloc(
                 file_size + 1
             );

             if (buffer == NULL)
             {
                 fclose(
                     file
                 );

                 return -1;
             }

             size_t actually_read =
             fread(
                 buffer,
                 1,
                 file_size,
                 file
             );

             fclose(
                 file
             );

             if (actually_read !=
                 file_size)
             {
                 free(
                     buffer
                 );

                 fprintf(
                     stderr,
                     "ERROR: cannot read configuration file '%s'.\n",
                     filename
                 );

                 return -1;
             }

             buffer[file_size] =
             '\0';

             ConfigArguments arguments;

             config_arguments_init(
                 &arguments
             );

             if (config_arguments_push(
                 &arguments,
                 program_name,
                 strlen(program_name)
             ) != 0)
             {
                 free(
                     buffer
                 );

                 config_arguments_free(
                     &arguments
                 );

                 return -1;
             }

             size_t i =
             0;

             while (i < file_size)
             {
                 /*
                  * Skip whitespace and full-line comments.
                  */
                 for (;;)
                 {
                     while (i < file_size &&
                         (
                             buffer[i] == ' ' ||
                             buffer[i] == '\t' ||
                             buffer[i] == '\r' ||
                             buffer[i] == '\n'
                         ))
                     {
                         ++i;
                     }

                     if (i < file_size &&
                         buffer[i] == '#')
                     {
                         while (i < file_size &&
                             buffer[i] != '\n')
                         {
                             ++i;
                         }

                         continue;
                     }

                     break;
                 }

                 if (i >= file_size)
                 {
                     break;
                 }

                 /*
                  * Build one token. Quotes are removed from the value.
                  */
                 char quote =
                 '\0';

                 if (buffer[i] == '"' ||
                     buffer[i] == '\'')
                 {
                     quote =
                     buffer[i];

                     ++i;
                 }

                 size_t token_capacity =
                 64;

                 size_t token_length =
                 0;

                 char *token =
                 malloc(
                     token_capacity
                 );

                 if (token == NULL)
                 {
                     free(
                         buffer
                     );

                     config_arguments_free(
                         &arguments
                     );

                     return -1;
                 }

                 while (i < file_size)
                 {
                     char c =
                     buffer[i];

                     if (quote != '\0')
                     {
                         if (c == quote)
                         {
                             ++i;
                             break;
                         }
                     }
                     else
                     {
                         if (c == ' ' ||
                             c == '\t' ||
                             c == '\r' ||
                             c == '\n')
                         {
                             break;
                         }

                         if (c == '#')
                         {
                             break;
                         }
                     }

                     if (token_length + 1 >=
                         token_capacity)
                     {
                         token_capacity *=
                         2;

                         char *new_token =
                         realloc(
                             token,
                             token_capacity
                         );

                         if (new_token == NULL)
                         {
                             free(
                                 token
                             );

                             free(
                                 buffer
                             );

                             config_arguments_free(
                                 &arguments
                             );

                             return -1;
                         }

                         token =
                         new_token;
                     }

                     token[
                         token_length
                     ] =
                     c;

                     ++token_length;
                     ++i;
                 }

                 if (quote != '\0' &&
                     (
                         i == file_size &&
                         (
                             file_size == 0 ||
                             buffer[file_size - 1] != quote
                         )
                     ))
                 {
                     fprintf(
                         stderr,
                         "ERROR: unterminated quote in configuration file '%s'.\n",
                         filename
                     );

                     free(
                         token
                     );

                     free(
                         buffer
                     );

                     config_arguments_free(
                         &arguments
                     );

                     return -1;
                 }

                 token[
                     token_length
                 ] =
                 '\0';

                 if (token_length > 0)
                 {
                     if (config_arguments_push(
                         &arguments,
                         token,
                         token_length
                     ) != 0)
                     {
                         free(
                             token
                         );

                         free(
                             buffer
                         );

                         config_arguments_free(
                             &arguments
                         );

                         return -1;
                     }
                 }

                 free(
                     token
                 );

                 /*
                  * Skip the separator after an unquoted token.
                  *
                  * If '#' begins a comment, skip to end of line.
                  */
                 if (quote == '\0' &&
                     i < file_size &&
                     buffer[i] == '#')
                 {
                     while (i < file_size &&
                         buffer[i] != '\n')
                     {
                         ++i;
                     }
                 }
             }

             free(
                 buffer
             );

             int result =
             config_parse_cli(
                 config,
                 arguments.argc,
                 arguments.argv
             );

             config_arguments_free(
                 &arguments
             );

             return result;
         }


         static int load_runtime_configuration(
             PeecConfig *config,
             int argc,
             char **argv
         )
         {
             if (config == NULL ||
                 argv == NULL ||
                 argc < 1)
             {
                 return -1;
             }

             if (argc > 1 && strcmp(argv[1], "--config") == 0) {
                 if (argc != 3) {
                     fprintf(stderr, "ERROR: usage: %s --config FILE\n", argv[0]);
                     return -1;
                 }
                 printf("Configuration file   : %s\n", argv[2]);
                 return config_parse_file_as_cli(config, argv[0], argv[2]);
             }

             if (argc > 1 || PEEC_INPUT_MODE ==
                 PEEC_INPUT_COMMAND_LINE)
             {
                 printf(
                     "Configuration source : COMMAND LINE\n"
                 );

                 return config_parse_cli(
                     config,
                     argc,
                     argv
                 );
             }

             if (PEEC_INPUT_MODE ==
                 PEEC_INPUT_FILE)
             {
                 printf(
                     "Configuration source : FILE\n"
                 );

                 printf(
                     "Configuration file   : %s\n",
                     PEEC_CONFIG_FILE
                 );

                 return config_parse_file_as_cli(
                     config,
                     argv[0],
                     PEEC_CONFIG_FILE
                 );
             }

             fprintf(
                 stderr,
                 "ERROR: unknown PEEC input mode.\n"
             );

             return -1;
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
             PeecConfig config;

             config_set_defaults(
                 &config
             );

             int parse_result =
             load_runtime_configuration(
                 &config,
                 argc,
                 argv
             );

             if (parse_result > 0)
             {
                 config_free(&config);
                 return EXIT_SUCCESS;
             }

             if (parse_result < 0)
             {
                 config_print_help(
                     argv[0]
                 );

                 config_free(&config);
                 return EXIT_FAILURE;
             }

             omp_set_dynamic(
                 0
             );

             omp_set_num_threads(
                 config.parallel_threads
             );

             printf("\n");
             printf("============================================================\n");
             printf("PEEC / PARTIAL EQUIVALENT CIRCUIT\n");
             printf("============================================================\n");

             printf(
                 "Task               : %s\n",
                 task_name(
                     config.task
                 )
             );

             printf(
                 "Physics            : %s\n",
                 physics_name(
                     config.physics
                 )
             );

             printf(
                 "Parallel           : %d OpenMP thread(s)\n",
                    config.parallel_threads
             );

             printf(
                 "Save matrices      : %s\n",
                 config.save_matrices
                 ? "YES"
                 : "NO"
             );

             printf(
                 "ParaView output    : %s\n",
                 config.write_vtk
                 ? "YES"
                 : "NO"
             );

             if (config.task == PEEC_TASK_LIGHTNING)
             {
                 printf(
                     "Time order         : %d (%s)\n",
                        config.time_order,
                        config.time_order == 1
                        ? "Backward Euler"
                        : "trapezoidal"
                 );

                 printf(
                     "Excitation         : %s\n",
                     excitation_name(
                         config.excitation
                     )
                 );
             }

             printf("============================================================\n");


             /*
              * Retarded mode is kept in config because it is part of the
              * final architecture, but the current L/P/harmonic/transient
              * kernels in this main are quasistatic.
              */
             if (config.physics ==
                 PEEC_PHYSICS_RETARDED)
             {
                 fprintf(
                     stderr,
                     "ERROR: --physics retarded is not connected to the "
                     "final driver yet.\n"
                 );

                 config_free(&config);
                 return EXIT_FAILURE;
             }


             int status =
             -1;

             switch (config.task)
             {
                 case PEEC_TASK_MESH_INFO:
                 {
                     status =
                     run_mesh_info(
                         &config
                     );

                     break;
                 }

                 case PEEC_TASK_SCATTERING:
                 {
                     status =
                     run_scattering(
                         &config
                     );

                     break;
                 }

                 case PEEC_TASK_RCS:
                 {
                     status =
                     run_rcs(
                         &config
                     );

                     break;
                 }

                 case PEEC_TASK_LIGHTNING:
                 {
                     status =
                     run_lightning(
                         &config
                     );

                     break;
                 }

                 case PEEC_TASK_NONE:
                 default:
                 {
                     fprintf(
                         stderr,
                         "ERROR: invalid task.\n"
                     );

                     status =
                     -1;

                     break;
                 }
             }

             if (status != 0)
             {
                 fprintf(
                     stderr,
                     "\nPEEC calculation failed.\n"
                 );

                 config_free(&config);
                 return EXIT_FAILURE;
             }

             printf("\nDone.\n");

             config_free(&config);
                 return EXIT_SUCCESS;
         }
