#ifndef PEEC_INCIDENT_FIELD_H
#define PEEC_INCIDENT_FIELD_H

#include "complex_matrix.h"


/*
 * ============================================================
 * ПОЛЯРИЗАЦИЯ ПАДАЮЩЕЙ ВОЛНЫ
 * ============================================================
 *
 * Сохраняем соглашение, аналогичное исходному Fortran-коду:
 *
 *     pol = 1 -> вертикальная поляризация e_theta
 *
 *     pol = 0 -> горизонтальная поляризация e_phi.
 */
typedef enum
{
    INCIDENT_POLARIZATION_HORIZONTAL = 0,
    INCIDENT_POLARIZATION_VERTICAL = 1

} IncidentPolarization;


/*
 * ============================================================
 * ТИП ПАДАЮЩЕГО ПОЛЯ
 * ============================================================
 */
typedef enum
{
    /*
     * Реальная монохроматическая волна:
     *
     *     E(r,t) =
     *
     *     E0 cos(
     *         omega*t
     *         - k0*k_hat.r
     *         + phase
     *     ).
     */
    INCIDENT_WAVE_COS = 0,

    /*
     * Реальная монохроматическая волна:
     *
     *     E(r,t) =
     *
     *     E0 sin(
     *         omega*t
     *         - k0*k_hat.r
     *         + phase
     *     ).
     */
    INCIDENT_WAVE_SIN = 1,

    /*
     * Комплексное представление.
     *
     * Используем соглашение:
     *
     *     exp(j*omega*t).
     *
     * Тогда пространственный фазор:
     *
     *     E_hat(r) =
     *
     *     E0 *
     *     exp(
     *         -j*k0*k_hat.r
     *         +j*phase
     *     ).
     *
     * Физическое поле:
     *
     *     E(r,t) =
     *
     *     Re{
     *         E_hat(r) exp(j*omega*t)
     *     }.
     */
    INCIDENT_WAVE_COMPLEX = 2,

    /*
     * Импульсная плоская волна:
     *
     *     E_inc(tau) =
     *
     *     k *
     *     (
     *         exp(-alpha*tau)
     *         -
     *         exp(-beta*tau)
     *     ).
     *
     * Здесь:
     *
     *     tau =
     *
     *     t - k_hat.r / c0.
     *
     * При tau < 0 поле равно нулю.
     */
    INCIDENT_WAVE_PULSE = 3

} IncidentWaveType;


/*
 * ============================================================
 * ПАРАМЕТРЫ ПАДАЮЩЕГО ПОЛЯ
 * ============================================================
 */
typedef struct
{
    /*
     * Тип волны.
     */
    IncidentWaveType wave_type;

    /*
     * Поляризация.
     */
    IncidentPolarization polarization;


    /*
     * --------------------------------------------------------
     * УГЛЫ
     * --------------------------------------------------------
     *
     * В градусах и радианах.
     */
    double theta_deg;
    double phi_deg;

    double theta;
    double phi;


    /*
     * --------------------------------------------------------
     * НАПРАВЛЕНИЕ РАСПРОСТРАНЕНИЯ
     * --------------------------------------------------------
     *
     * Согласно принятому соглашению:
     *
     *     k_hat =
     *
     *     -(
     *         sin(theta) cos(phi),
     *         sin(theta) sin(phi),
     *         cos(theta)
     *     ).
     *
     * Это единичный вектор.
     */
    double k_hat[3];


    /*
     * --------------------------------------------------------
     * БАЗИС ПОЛЯРИЗАЦИИ
     * --------------------------------------------------------
     */
    double e_theta[3];
    double e_phi[3];


    /*
     * --------------------------------------------------------
     * АМПЛИТУДА ВЕКТОРА E0
     * --------------------------------------------------------
     *
     * Для гармонической волны:
     *
     *     E_amp [V/m].
     */
    double E_amp;

    /*
     * Вектор амплитуды электрического поля:
     *
     *     E0 = E_amp * polarization_vector.
     *
     * Единицы:
     *
     *     V/m.
     */
    double E0[3];


    /*
     * --------------------------------------------------------
     * ГАРМОНИЧЕСКИЕ ПАРАМЕТРЫ
     * --------------------------------------------------------
     *
     * frequency [Hz].
     */
    double frequency;

    /*
     * omega = 2*pi*f [rad/s].
     */
    double omega;

    /*
     * k0 = omega/c0 [rad/m].
     *
     * Не путать с коэффициентом k импульса.
     */
    double wave_number;

    /*
     * Начальная фаза [rad].
     */
    double phase;


    /*
     * --------------------------------------------------------
     * ИМПУЛЬС
     * --------------------------------------------------------
     *
     *     E_inc(t) =
     *
     *     pulse_k *
     *     (
     *         exp(-alpha*t)
     *         -
     *         exp(-beta*t)
     *     ).
     *
     *
     * pulse_k имеет размерность:
     *
     *     V/m.
     *
     * Поэтому для импульсного режима E_amp отдельно
     * не используется как амплитуда.
     */
    double pulse_k;

    /*
     * alpha [1/s].
     */
    double pulse_alpha;

    /*
     * beta [1/s].
     */
    double pulse_beta;

} IncidentField;


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void incident_field_init(
    IncidentField *field
);


/*
 * ============================================================
 * НАСТРОЙКА МОНОХРОМАТИЧЕСКОЙ ВОЛНЫ
 * ============================================================
 *
 * wave_type должен быть:
 *
 *     INCIDENT_WAVE_COS
 *     INCIDENT_WAVE_SIN
 *     INCIDENT_WAVE_COMPLEX.
 */
int incident_field_configure_harmonic(
    IncidentField *field,
    IncidentWaveType wave_type,
    IncidentPolarization polarization,
    double E_amp,
    double frequency,
    double theta_deg,
    double phi_deg,
    double phase
);


/*
 * ============================================================
 * НАСТРОЙКА ИМПУЛЬСНОЙ ВОЛНЫ
 * ============================================================
 *
 *     E_inc(t) =
 *
 *     pulse_k *
 *     (
 *         exp(-alpha*t)
 *         -
 *         exp(-beta*t)
 *     ).
 */
int incident_field_configure_pulse(
    IncidentField *field,
    IncidentPolarization polarization,
    double pulse_k,
    double alpha,
    double beta,
    double theta_deg,
    double phi_deg
);


/*
 * ============================================================
 * РЕАЛЬНОЕ ПОЛЕ E(r,t)
 * ============================================================
 *
 * Используется для:
 *
 *     COS
 *     SIN
 *     PULSE.
 *
 * Результат:
 *
 *     E[3] [V/m].
 */
int incident_field_evaluate_time(
    const IncidentField *field,
    const double r[3],
    double time,
    double E[3]
);


/*
 * ============================================================
 * КОМПЛЕКСНЫЙ ФАЗОР E_hat(r)
 * ============================================================
 *
 * Используется для частотной PEEC-задачи.
 *
 * Соглашение:
 *
 *     E(r,t) =
 *
 *     Re{
 *         E_hat(r) exp(j*omega*t)
 *     }.
 *
 *
 * Для COS:
 *
 *     E_hat =
 *
 *     E0 exp(
 *         -j*k0*k_hat.r
 *         +j*phase
 *     ).
 *
 *
 * Для SIN автоматически учитывается:
 *
 *     sin(x) = Re{-j exp(j*x)}.
 *
 *
 * Для COMPLEX используется непосредственно
 * комплексная экспонента.
 */
int incident_field_evaluate_phasor(
    const IncidentField *field,
    const double r[3],
    Complex E[3]
);


/*
 * Диагностика.
 */
void incident_field_print_info(
    const IncidentField *field
);


#endif