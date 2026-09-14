#ifndef PEEC_CONFIG_H
#define PEEC_CONFIG_H

/*
 * ============================================================
 * ТИП РЕШАЕМОЙ ЗАДАЧИ
 * ============================================================
 */
typedef enum
{
    PEEC_TASK_NONE = 0,
    PEEC_TASK_MESH_INFO,
    PEEC_TASK_SCATTERING,
    PEEC_TASK_RCS,
    PEEC_TASK_LIGHTNING

} PeecTask;


/*
 * ============================================================
 * ФИЗИЧЕСКАЯ МОДЕЛЬ
 * ============================================================
 */
typedef enum
{
    PEEC_PHYSICS_QUASISTATIC = 0,
    PEEC_PHYSICS_RETARDED

} PeecPhysicsMode;


/*
 * ============================================================
 * ПОЛЯРИЗАЦИЯ ПАДАЮЩЕЙ ВОЛНЫ
 * ============================================================
 */
typedef enum
{
    PEEC_POLARIZATION_VERTICAL = 0,
    PEEC_POLARIZATION_HORIZONTAL

} PeecPolarization;


/*
 * ============================================================
 * ТИП LIGHTNING-ЗАДАЧИ
 * ============================================================
 *
 * CLOSED:
 *
 *     один полностью закрытый корпус;
 *
 *     выполняется только transient Stage 1.
 *
 *
 * TWO_STAGE:
 *
 *     Stage 1:
 *
 *         закрытая PEC-модель корпуса;
 *
 *     Stage 2:
 *
 *         открытая модель корпуса с апертурой /
 *         внутренней неоднородностью.
 */
typedef enum
{
    PEEC_LIGHTNING_CLOSED_BODY = 0,
    PEEC_LIGHTNING_TWO_STAGE,
    PEEC_LIGHTNING_SLOT_CELLS

} PeecLightningCase;


/*
 * ============================================================
 * TRANSIENT EXCITATION
 * ============================================================
 *
 * LIGHTNING_CURRENT:
 *
 *     узловой ток молнии между strike / return.
 *
 *
 * INCIDENT_PULSE:
 *
 *     внешняя падающая импульсная плоская волна:
 *
 *         E(tau) =
 *
 *         K [
 *             exp(-alpha*tau)
 *             -
 *             exp(-beta*tau)
 *           ],
 *
 *     tau = t - k_hat.r/c0.
 */
typedef enum
{
    PEEC_EXCITATION_LIGHTNING_CURRENT = 0,
    PEEC_EXCITATION_INCIDENT_PULSE

} PeecExcitation;


/*
 * ============================================================
 * ОБЩАЯ КОНФИГУРАЦИЯ
 * ============================================================
 */
typedef struct
{
    /*
     * --------------------------------------------------------
     * ОБЩИЕ ПАРАМЕТРЫ
     * --------------------------------------------------------
     */
    PeecTask task;
    PeecPhysicsMode physics;

    /*
     * Основная сетка.
     *
     * Используется для:
     *
     *     mesh-info
     *     scattering
     *     rcs.
     *
     *
     * Для lightning может использоваться как сокращение
     * для closed_mesh_file в режиме closed.
     */
    const char *mesh_file;

    int parallel_threads;


    /*
     * --------------------------------------------------------
     * SCATTERING / RCS
     * --------------------------------------------------------
     */
    double frequency;
    double field_amplitude;

    double theta_deg;
    double phi_deg;
    double phase_rad;

    PeecPolarization polarization;


    /*
     * --------------------------------------------------------
     * RCS
     * --------------------------------------------------------
     *
     * ka:
     *
     *     безразмерный электрический размер.
     *
     *
     * characteristic_length:
     *
     *     характерный размер a [m].
     *
     *
     * frequency для RCS может быть вычислена как:
     *
     *     k = ka / a
     *
     *     f = c k / (2*pi).
     */
    double ka;
    double characteristic_length;

    int rcs_use_dual;
    int rcs_order;

    double rcs_observation_theta;
    double rcs_phi_start;
    double rcs_phi_end;
    double rcs_phi_step;

    const char *output_file;


    /*
     * --------------------------------------------------------
     * OUTPUT / CACHE
     * --------------------------------------------------------
     *
     * Matrix Market:
     *
     *     L.mtx
     *     A.mtx
     *     P.mtx
     *     PAT.mtx
     *     R.mtx
     */
    int save_matrices;
    const char *matrix_directory;

    /*
     * ParaView Legacy VTK.
     */
    int write_vtk;
    const char *vtk_directory;

    /*
     * Для transient:
     *
     * сохранять каждый vtk_every шаг.
     */
    int vtk_every;


    /*
     * --------------------------------------------------------
     * LIGHTNING
     * --------------------------------------------------------
     */
    PeecLightningCase lightning_case;

    /*
     * Тип transient-возбуждения.
     */
    PeecExcitation excitation;

    /*
     * Stage 1 / closed body mesh.
     */
    const char *closed_mesh_file;

    /*
     * Stage 2 / open body mesh.
     *
     * Нужна только для:
     *
     *     PEEC_LIGHTNING_TWO_STAGE.
     */
    const char *open_mesh_file;

    /*
     * Карта переноса Stage 1 -> Stage 2.
     *
     * Нужна только для двухэтапной задачи.
     */
    const char *aperture_map_file;


    /*
     * --------------------------------------------------------
     * ONE-STAGE DISTRIBUTED SLOT MODEL
     * --------------------------------------------------------
     *
     * Used only for:
     *
     *     PEEC_LIGHTNING_SLOT_CELLS.
     *
     * The open mesh is solved once and the aperture is replaced
     * by distributed equivalent LC cells.
     */
    const char *slot_cells_file;

    double slot_width;
    double slot_wall_span;
    double slot_wall_thickness;

    double slot_epsilon_r;
    double slot_mu_r;


    /*
     * Временная сетка.
     *
     * [s]
     */
    double dt;
    double t_end;

    /*
     * Порядок неявной аппроксимации по времени:
     *
     *     1 -> Backward Euler;
     *     2 -> trapezoidal rule.
     *
     * По умолчанию:
     *
     *     2.
     */
    int time_order;


    /*
     * Параметры тока молнии:
     *
     *     I_L(t)
     *
     *     =
     *
     *     K [
     *         exp(-alpha*(t-delay))
     *         -
     *         exp(-beta*(t-delay))
     *       ].
     */
    double lightning_K;
    double lightning_alpha;
    double lightning_beta;
    double lightning_delay;


    /*
     * --------------------------------------------------------
     * INCIDENT PULSE
     * --------------------------------------------------------
     *
     * Амплитуда pulse_K имеет размерность [V/m].
     *
     * pulse_alpha, pulse_beta:
     *
     *     [1/s].
     *
     * Направление и поляризация берутся из общих:
     *
     *     theta_deg
     *     phi_deg
     *     polarization.
     */
    double pulse_K;
    double pulse_alpha;
    double pulse_beta;


    /*
     * Координаты точек подключения молнии.
     *
     * [m]
     */
    double strike_xyz[3];
    double return_xyz[3];

    int strike_xyz_set;
    int return_xyz_set;

    /*
     * Прямое задание внутренних индексов узлов solver-а:
     *
     *     0 <= node < mesh.n_nodes.
     *
     * Это НЕ Gmsh node tag, а индекс массива mesh.xyz.
     */
    int strike_node;
    int return_node;

    int strike_node_set;
    int return_node_set;


    /*
     * --------------------------------------------------------
     * STAGE-2 INTERNAL SHUNT
     * --------------------------------------------------------
     */
    int use_internal_shunt;

    double shunt_resistance;

    double shunt_a_xyz[3];
    double shunt_b_xyz[3];

    int shunt_a_xyz_set;
    int shunt_b_xyz_set;

    /*
     * Прямое задание внутренних индексов узлов Stage 2.
     */
    int shunt_a_node;
    int shunt_b_node;

    int shunt_a_node_set;
    int shunt_b_node_set;

    /* Owned copies of option strings; release with config_free(). */
    struct PeecConfigString *owned_strings;
} PeecConfig;


/*
 * ============================================================
 * API
 * ============================================================
 */
/* Initialize once, parse as needed, then free on every exit path.
 * PeecConfig owns option strings and must not be freed via a shallow copy. */
void config_free(PeecConfig *config);

void config_set_defaults(
    PeecConfig *config
);

int config_parse_cli(
    PeecConfig *config,
    int argc,
    char **argv
);

void config_print_help(
    const char *program_name
);

#endif
