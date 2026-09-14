#include "resistance.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * ============================================================
 * ПАРАМЕТРЫ ПО УМОЛЧАНИЮ
 * ============================================================
 *
 * В электромагнитной задаче безопаснее по умолчанию считать
 * поверхность идеальным проводником, пока пользователь явно
 * не задаст конечное поверхностное сопротивление.
 */
void resistance_options_default(
    ResistanceOptions *options,
    int parallel_threads
)
{
    if (options == NULL)
    {
        return;
    }

    options->model =
    PEEC_CONDUCTOR_PEC;

    options->sheet_resistance =
    0.0;

    options->parallel_threads =
    parallel_threads;
}


/*
 * ============================================================
 * ИНИЦИАЛИЗАЦИЯ
 * ============================================================
 */
void resistance_init(
    ResistanceMatrix *R
)
{
    if (R == NULL)
    {
        return;
    }

    memset(
        R,
        0,
        sizeof(*R)
    );
}


/*
 * ============================================================
 * ОСВОБОЖДЕНИЕ ПАМЯТИ
 * ============================================================
 */
void resistance_free(
    ResistanceMatrix *R
)
{
    if (R == NULL)
    {
        return;
    }

    free(R->diagonal);

    resistance_init(R);
}


/*
 * ============================================================
 * РАСЧЕТ R
 * ============================================================
 */
int resistance_compute_matrix(
    const Mesh *mesh,
    const DualMesh *dual,
    const ResistanceOptions *options,
    ResistanceMatrix *R
)
{
    if (mesh == NULL ||
        dual == NULL ||
        options == NULL ||
        R == NULL ||
        options->parallel_threads < 1)
    {
        return -1;
    }

    if (mesh->n_edges == 0 ||
        dual->n_edges != mesh->n_edges ||
        mesh->edge_lengths == NULL ||
        dual->edge_region_areas == NULL)
    {
        fprintf(
            stderr,
            "ERROR: invalid mesh or dual mesh in resistance calculation.\n"
        );

        return -1;
    }


    /*
     * Для модели с конечным поверхностным сопротивлением:
     *
     *     R_s >= 0.
     */
    if (options->model == PEEC_CONDUCTOR_SHEET_RESISTANCE)
    {
        if (!isfinite(options->sheet_resistance) ||
            options->sheet_resistance < 0.0)
        {
            fprintf(
                stderr,
                "ERROR: invalid sheet resistance %.9e Ohm.\n",
                options->sheet_resistance
            );

            return -1;
        }
    }
    else if (options->model != PEEC_CONDUCTOR_PEC)
    {
        fprintf(
            stderr,
            "ERROR: unknown conductor model.\n"
        );

        return -1;
    }


    resistance_free(R);

    R->n =
    mesh->n_edges;

    R->diagonal =
    calloc(
        R->n,
        sizeof(double)
    );

    if (R->diagonal == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate resistance vector.\n"
        );

        resistance_init(R);

        return -1;
    }


    /*
     * ========================================================
     * PEC
     * ========================================================
     *
     * calloc уже создал:
     *
     *     R_i = 0.
     *
     * Поэтому дополнительных вычислений не требуется.
     */
    if (options->model == PEEC_CONDUCTOR_PEC)
    {
        return 0;
    }


    /*
     * ========================================================
     * ПОВЕРХНОСТНОЕ СОПРОТИВЛЕНИЕ
     * ========================================================
     *
     * Для каждого ребра:
     *
     *     R_i = R_s * l_i / w_i,
     *
     * где:
     *
     *     w_i = S_i^e / l_i.
     *
     * Следовательно:
     *
     *     R_i = R_s * l_i^2 / S_i^e.
     */
    double sheet_resistance =
    options->sheet_resistance;

    #pragma omp parallel for schedule(static) num_threads(options->parallel_threads)
    for (long long i_signed = 0;
         i_signed < (long long)mesh->n_edges;
    ++i_signed)
         {
             size_t i =
             (size_t)i_signed;

             double length =
             mesh->edge_lengths[i];

             double region_area =
             dual->edge_region_areas[i];


             /*
              * Геометрические величины должны быть положительными.
              *
              * Ошибки геометрии уже желательно ловить раньше,
              * однако здесь дополнительно защищаем вычисление.
              */
             if (length <= 0.0 ||
                 region_area <= 0.0 ||
                 !isfinite(length) ||
                 !isfinite(region_area))
             {
                 R->diagonal[i] =
                 NAN;

                 continue;
             }


             /*
              * Эффективная ширина поверхностного тока:
              *
              *     w_i = S_i^e / l_i.
              */
             double width =
             region_area / length;


             /*
              * Сопротивление реберной ветви:
              *
              *     R_i = R_s * l_i / w_i.
              */
             R->diagonal[i] =
             sheet_resistance
             * length
             / width;
         }


         /*
          * Проверяем результат после параллельного цикла.
          */
         for (size_t i = 0;
              i < R->n;
    ++i)
              {
                  if (!isfinite(R->diagonal[i]) ||
                      R->diagonal[i] < 0.0)
                  {
                      fprintf(
                          stderr,
                          "ERROR: invalid resistance at edge %zu: %.9e Ohm.\n",
                          i,
                          R->diagonal[i]
                      );

                      resistance_free(R);

                      return -1;
                  }
              }

              return 0;
}


/*
 * ============================================================
 * ДОСТУП К R_i
 * ============================================================
 */
double resistance_get(
    const ResistanceMatrix *R,
    size_t edge
)
{
    if (R == NULL ||
        R->diagonal == NULL ||
        edge >= R->n)
    {
        return NAN;
    }

    return R->diagonal[edge];
}


/*
 * ============================================================
 * ДИАГНОСТИКА
 * ============================================================
 */
void resistance_print_info(
    const ResistanceMatrix *R,
    const ResistanceOptions *options
)
{
    if (R == NULL ||
        options == NULL)
    {
        return;
    }

    printf("\n");
    printf("========================================\n");
    printf("RESISTANCE MATRIX\n");
    printf("========================================\n");

    printf(
        "Edges              : %zu\n",
        R->n
    );


    if (options->model == PEEC_CONDUCTOR_PEC)
    {
        printf(
            "Conductor model    : PEC\n"
        );

        printf(
            "Sheet resistance   : 0 Ohm\n"
        );
    }
    else if (options->model == PEEC_CONDUCTOR_SHEET_RESISTANCE)
    {
        printf(
            "Conductor model    : sheet resistance\n"
        );

        printf(
            "Sheet resistance   : %.9e Ohm\n",
            options->sheet_resistance
        );
    }
    else
    {
        printf(
            "Conductor model    : unknown\n"
        );
    }


    if (R->diagonal == NULL ||
        R->n == 0)
    {
        printf(
            "Matrix             : not allocated\n"
        );

        printf(
            "========================================\n"
        );

        return;
    }


    double min_value =
    R->diagonal[0];

    double max_value =
    R->diagonal[0];

    for (size_t i = 1;i < R->n;++i)
         {
             double value =
             R->diagonal[i];

             if (value < min_value)
             {
                 min_value =
                 value;
             }

             if (value > max_value)
             {
                 max_value =
                 value;
             }
         }


         double memory_kib =
         ((double)R->n
         * sizeof(double))
         / 1024.0;


         printf(
             "Minimum R_i        : %.9e Ohm\n",
             min_value
         );

         printf(
             "Maximum R_i        : %.9e Ohm\n",
             max_value
         );

         printf(
             "Diagonal memory    : %.3f KiB\n",
             memory_kib
         );

         printf(
             "========================================\n"
         );
}
