#include "incident_field.h"

#include <math.h>
#include <stdio.h>
#include <string.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


/*
 * Скорость света в вакууме:
 *
 *     c0 [m/s].
 */
static const double C0 =
    299792458.0;


/*
 * ============================================================
 * ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
 * ============================================================
 */

static double vector_dot(
    const double a[3],
    const double b[3]
)
{
    return
        a[0] * b[0]
        + a[1] * b[1]
        + a[2] * b[2];
}


static double vector_norm(
    const double a[3]
)
{
    return sqrt(
        a[0] * a[0]
        + a[1] * a[1]
        + a[2] * a[2]
    );
}


static int vector_normalize(
    double a[3]
)
{
    double nrm =
        vector_norm(a);

    if (!isfinite(nrm) ||
        nrm <= 1.0e-14)
    {
        return -1;
    }

    a[0] /= nrm;
    a[1] /= nrm;
    a[2] /= nrm;

    return 0;
}


/*
 * ============================================================
 * ПОСТРОЕНИЕ НАПРАВЛЕНИЯ И ПОЛЯРИЗАЦИИ
 * ============================================================
 *
 * Полностью повторяем принятое пользователем соглашение:
 *
 *     k_hat =
 *
 *     -(
 *         sin(theta) cos(phi),
 *         sin(theta) sin(phi),
 *         cos(theta)
 *     ).
 *
 *
 * Вертикальная поляризация:
 *
 *     e_theta =
 *
 *     (
 *         cos(theta) cos(phi),
 *         cos(theta) sin(phi),
 *         -sin(theta)
 *     ).
 *
 *
 * Горизонтальная:
 *
 *     e_phi =
 *
 *     (
 *         -sin(phi),
 *          cos(phi),
 *          0
 *     ).
 */
static int incident_field_build_basis(
    IncidentField *field
)
{
    if (field == NULL)
    {
        return -1;
    }

    double theta =
        field->theta;

    double phi =
        field->phi;


    /*
     * ========================================================
     * НАПРАВЛЕНИЕ РАСПРОСТРАНЕНИЯ
     * ========================================================
     */
    field->k_hat[0] =
        -sin(theta) * cos(phi);

    field->k_hat[1] =
        -sin(theta) * sin(phi);

    field->k_hat[2] =
        -cos(theta);


    /*
     * ========================================================
     * ВЕРТИКАЛЬНЫЙ БАЗИС
     * ========================================================
     */
    field->e_theta[0] =
        cos(theta) * cos(phi);

    field->e_theta[1] =
        cos(theta) * sin(phi);

    field->e_theta[2] =
        -sin(theta);


    /*
     * ========================================================
     * ГОРИЗОНТАЛЬНЫЙ БАЗИС
     * ========================================================
     */
    field->e_phi[0] =
        -sin(phi);

    field->e_phi[1] =
        cos(phi);

    field->e_phi[2] =
        0.0;


    /*
     * На всякий случай нормируем,
     * как это делалось в исходном Fortran-коде.
     */
    if (vector_normalize(field->k_hat) != 0 ||
        vector_normalize(field->e_theta) != 0 ||
        vector_normalize(field->e_phi) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot build incident-field basis.\n"
        );

        return -1;
    }


    /*
     * Проверяем поперечность:
     *
     *     k_hat . e_theta = 0
     *     k_hat . e_phi   = 0.
     */
    double dot_theta =
        fabs(
            vector_dot(
                field->k_hat,
                field->e_theta
            )
        );

    double dot_phi =
        fabs(
            vector_dot(
                field->k_hat,
                field->e_phi
            )
        );

    if (dot_theta > 1.0e-12 ||
        dot_phi > 1.0e-12)
    {
        fprintf(
            stderr,
            "ERROR: incident polarization is not transverse.\n"
        );

        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ПОСТРОЕНИЕ E0
 * ============================================================
 */
static int incident_field_build_E0(
    IncidentField *field
)
{
    if (field == NULL)
    {
        return -1;
    }

    const double *polarization_vector;

    if (field->polarization ==
        INCIDENT_POLARIZATION_VERTICAL)
    {
        polarization_vector =
            field->e_theta;
    }
    else if (field->polarization ==
             INCIDENT_POLARIZATION_HORIZONTAL)
    {
        polarization_vector =
            field->e_phi;
    }
    else
    {
        fprintf(
            stderr,
            "ERROR: unknown incident polarization.\n"
        );

        return -1;
    }

    field->E0[0] =
        field->E_amp
        * polarization_vector[0];

    field->E0[1] =
        field->E_amp
        * polarization_vector[1];

    field->E0[2] =
        field->E_amp
        * polarization_vector[2];

    return 0;
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void incident_field_init(
    IncidentField *field
)
{
    if (field == NULL)
    {
        return;
    }

    memset(
        field,
        0,
        sizeof(*field)
    );

    field->wave_type =
        INCIDENT_WAVE_COS;

    field->polarization =
        INCIDENT_POLARIZATION_VERTICAL;
}


/*
 * ============================================================
 * МОНОХРОМАТИЧЕСКАЯ ВОЛНА
 * ============================================================
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
)
{
    if (field == NULL)
    {
        return -1;
    }

    if (wave_type != INCIDENT_WAVE_COS &&
        wave_type != INCIDENT_WAVE_SIN &&
        wave_type != INCIDENT_WAVE_COMPLEX)
    {
        fprintf(
            stderr,
            "ERROR: invalid harmonic incident wave type.\n"
        );

        return -1;
    }

    if (!isfinite(E_amp) ||
        !isfinite(frequency) ||
        frequency <= 0.0 ||
        !isfinite(theta_deg) ||
        !isfinite(phi_deg) ||
        !isfinite(phase))
    {
        fprintf(
            stderr,
            "ERROR: invalid harmonic incident-field parameters.\n"
        );

        return -1;
    }


    incident_field_init(field);

    field->wave_type =
        wave_type;

    field->polarization =
        polarization;

    field->E_amp =
        E_amp;

    field->frequency =
        frequency;

    field->theta_deg =
        theta_deg;

    field->phi_deg =
        phi_deg;

    field->phase =
        phase;


    field->theta =
        theta_deg
        * M_PI
        / 180.0;

    field->phi =
        phi_deg
        * M_PI
        / 180.0;


    field->omega =
        2.0
        * M_PI
        * frequency;

    field->wave_number =
        field->omega
        / C0;


    if (incident_field_build_basis(
            field
        ) != 0)
    {
        return -1;
    }


    if (incident_field_build_E0(
            field
        ) != 0)
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ИМПУЛЬС
 * ============================================================
 */
int incident_field_configure_pulse(
    IncidentField *field,
    IncidentPolarization polarization,
    double pulse_k,
    double alpha,
    double beta,
    double theta_deg,
    double phi_deg
)
{
    if (field == NULL)
    {
        return -1;
    }

    if (!isfinite(pulse_k) ||
        !isfinite(alpha) ||
        !isfinite(beta) ||
        alpha <= 0.0 ||
        beta <= 0.0 ||
        !isfinite(theta_deg) ||
        !isfinite(phi_deg))
    {
        fprintf(
            stderr,
            "ERROR: invalid pulse incident-field parameters.\n"
        );

        return -1;
    }


    incident_field_init(field);

    field->wave_type =
        INCIDENT_WAVE_PULSE;

    field->polarization =
        polarization;

    field->pulse_k =
        pulse_k;

    field->pulse_alpha =
        alpha;

    field->pulse_beta =
        beta;

    field->theta_deg =
        theta_deg;

    field->phi_deg =
        phi_deg;


    field->theta =
        theta_deg
        * M_PI
        / 180.0;

    field->phi =
        phi_deg
        * M_PI
        / 180.0;


    if (incident_field_build_basis(
            field
        ) != 0)
    {
        return -1;
    }


    /*
     * Для импульса E0 используется как единичный
     * вектор поляризации.
     *
     * Амплитуда задается через pulse_k.
     */
    field->E_amp =
        1.0;

    if (incident_field_build_E0(
            field
        ) != 0)
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * ВРЕМЕННОЕ ПОЛЕ
 * ============================================================
 */
int incident_field_evaluate_time(
    const IncidentField *field,
    const double r[3],
    double time,
    double E[3]
)
{
    if (field == NULL ||
        r == NULL ||
        E == NULL ||
        !isfinite(time))
    {
        return -1;
    }


    double kr =
        vector_dot(
            field->k_hat,
            r
        );


    /*
     * ========================================================
     * COS / SIN
     * ========================================================
     */
    if (field->wave_type == INCIDENT_WAVE_COS ||
        field->wave_type == INCIDENT_WAVE_SIN)
    {
        /*
         * Фаза:
         *
         *     psi =
         *
         *     omega*t
         *     - k0*k_hat.r
         *     + phase.
         *
         *
         * Это согласовано с:
         *
         *     exp(j*omega*t)
         *
         * в частотной области.
         */
        double psi =
            field->omega * time
            - field->wave_number * kr
            + field->phase;

        double value;

        if (field->wave_type ==
            INCIDENT_WAVE_COS)
        {
            value =
                cos(psi);
        }
        else
        {
            value =
                sin(psi);
        }

        E[0] =
            field->E0[0] * value;

        E[1] =
            field->E0[1] * value;

        E[2] =
            field->E0[2] * value;

        return 0;
    }


    /*
     * ========================================================
     * ИМПУЛЬСНАЯ ВОЛНА
     * ========================================================
     *
     * В точке r импульс приходит с задержкой:
     *
     *     tau =
     *
     *     t - k_hat.r/c0.
     *
     *
     * Огибающая:
     *
     *     E_inc(tau) =
     *
     *     k *
     *     (
     *         exp(-alpha*tau)
     *         -
     *         exp(-beta*tau)
     *     ).
     */
    if (field->wave_type ==
        INCIDENT_WAVE_PULSE)
    {
        double tau =
            time
            - kr / C0;

        if (tau < 0.0)
        {
            E[0] = 0.0;
            E[1] = 0.0;
            E[2] = 0.0;

            return 0;
        }

        double envelope =
            field->pulse_k
            * (
                exp(
                    -field->pulse_alpha
                    * tau
                )
                -
                exp(
                    -field->pulse_beta
                    * tau
                )
            );

        E[0] =
            field->E0[0]
            * envelope;

        E[1] =
            field->E0[1]
            * envelope;

        E[2] =
            field->E0[2]
            * envelope;

        return 0;
    }


    /*
     * Для COMPLEX временное поле при необходимости
     * можно получать через phasor.
     *
     * Но основная функция для него:
     *
     *     incident_field_evaluate_phasor().
     */
    if (field->wave_type ==
        INCIDENT_WAVE_COMPLEX)
    {
        Complex phasor[3];

        if (incident_field_evaluate_phasor(
                field,
                r,
                phasor
            ) != 0)
        {
            return -1;
        }

        /*
         * Физическое поле:
         *
         *     Re{
         *         E_hat exp(j*omega*t)
         *     }.
         */
        double c =
            cos(
                field->omega
                * time
            );

        double s =
            sin(
                field->omega
                * time
            );

        for (int k = 0; k < 3; ++k)
        {
            E[k] =
                phasor[k].re * c
                - phasor[k].im * s;
        }

        return 0;
    }


    return -1;
}


/*
 * ============================================================
 * КОМПЛЕКСНОЕ ПОЛЕ
 * ============================================================
 *
 * Соглашение:
 *
 *     E(r,t) =
 *
 *     Re{
 *         E_hat(r)
 *         exp(j*omega*t)
 *     }.
 *
 *
 * Поэтому:
 *
 *     E_hat(r) =
 *
 *     E0 *
 *     exp(
 *         -j*k0*k_hat.r
 *         +j*phase
 *     ).
 */
int incident_field_evaluate_phasor(
    const IncidentField *field,
    const double r[3],
    Complex E[3]
)
{
    if (field == NULL ||
        r == NULL ||
        E == NULL)
    {
        return -1;
    }

    if (field->wave_type != INCIDENT_WAVE_COS &&
        field->wave_type != INCIDENT_WAVE_SIN &&
        field->wave_type != INCIDENT_WAVE_COMPLEX)
    {
        fprintf(
            stderr,
            "ERROR: phasor requested for non-harmonic field.\n"
        );

        return -1;
    }


    double kr =
        vector_dot(
            field->k_hat,
            r
        );


    /*
     * Пространственная фаза:
     *
     *     gamma =
     *
     *     -k0*k_hat.r
     *     + phase.
     */
    double gamma =
        -field->wave_number
        * kr
        + field->phase;


    /*
     * Для COS:
     *
     *     cos(psi)
     *
     *     = Re{
     *         exp(j*psi)
     *       }.
     *
     *
     * Для SIN:
     *
     *     sin(psi)
     *
     *     = Re{
     *         -j exp(j*psi)
     *       }.
     *
     * Поэтому SIN получает дополнительный:
     *
     *     -pi/2
     *
     * в фазоре.
     */
    if (field->wave_type ==
        INCIDENT_WAVE_SIN)
    {
        gamma -=
            M_PI / 2.0;
    }


    double c =
        cos(gamma);

    double s =
        sin(gamma);


    for (int k = 0; k < 3; ++k)
    {
        E[k].re =
            field->E0[k]
            * c;

        E[k].im =
            field->E0[k]
            * s;
    }

    return 0;
}


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void incident_field_print_info(
    const IncidentField *field
)
{
    if (field == NULL)
    {
        return;
    }

    printf("\n");
    printf("============================================================\n");
    printf("INCIDENT FIELD\n");
    printf("============================================================\n");

    printf(
        "theta       : %.6f deg\n",
        field->theta_deg
    );

    printf(
        "phi         : %.6f deg\n",
        field->phi_deg
    );

    printf(
        "k_hat       : % .9e  % .9e  % .9e\n",
        field->k_hat[0],
        field->k_hat[1],
        field->k_hat[2]
    );


    if (field->polarization ==
        INCIDENT_POLARIZATION_VERTICAL)
    {
        printf(
            "Polarization: VERTICAL (e_theta)\n"
        );
    }
    else
    {
        printf(
            "Polarization: HORIZONTAL (e_phi)\n"
        );
    }


    printf(
        "e_theta     : % .9e  % .9e  % .9e\n",
        field->e_theta[0],
        field->e_theta[1],
        field->e_theta[2]
    );

    printf(
        "e_phi       : % .9e  % .9e  % .9e\n",
        field->e_phi[0],
        field->e_phi[1],
        field->e_phi[2]
    );


    if (field->wave_type ==
        INCIDENT_WAVE_PULSE)
    {
        printf(
            "Wave type    : PULSE\n"
        );

        printf(
            "pulse k      : %.9e V/m\n",
            field->pulse_k
        );

        printf(
            "alpha        : %.9e 1/s\n",
            field->pulse_alpha
        );

        printf(
            "beta         : %.9e 1/s\n",
            field->pulse_beta
        );

        printf(
            "E direction  : % .9e  % .9e  % .9e\n",
            field->E0[0],
            field->E0[1],
            field->E0[2]
        );
    }
    else
    {
        printf(
            "E amplitude  : %.9e V/m\n",
            field->E_amp
        );

        printf(
            "E0           : % .9e  % .9e  % .9e V/m\n",
            field->E0[0],
            field->E0[1],
            field->E0[2]
        );

        printf(
            "Frequency    : %.9e Hz\n",
            field->frequency
        );

        printf(
            "omega        : %.9e rad/s\n",
            field->omega
        );

        printf(
            "wave number  : %.9e rad/m\n",
            field->wave_number
        );

        printf(
            "phase        : %.9e rad\n",
            field->phase
        );


        if (field->wave_type ==
            INCIDENT_WAVE_COS)
        {
            printf(
                "Wave type    : COS\n"
            );
        }
        else if (field->wave_type ==
                 INCIDENT_WAVE_SIN)
        {
            printf(
                "Wave type    : SIN\n"
            );
        }
        else
        {
            printf(
                "Wave type    : COMPLEX exp(j*omega*t)\n"
            );
        }
    }

    printf("============================================================\n");
}