#include "incidence.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void incidence_init(IncidenceMatrix *A)
{
    if (A == NULL)
    {
        return;
    }

    memset(A, 0, sizeof(*A));
}


void incidence_free(IncidenceMatrix *A)
{
    if (A == NULL)
    {
        return;
    }

    free(A->data);

    incidence_init(A);
}


int incidence_build(const Mesh *mesh, IncidenceMatrix *A)
{
    if (mesh == NULL || A == NULL)
    {
        return -1;
    }

    if (mesh->edges == NULL || mesh->n_edges == 0 || mesh->n_nodes == 0)
    {
        fprintf(stderr, "ERROR: mesh topology is not available.\n");
        return -1;
    }

    incidence_init(A);

    /*
     * Размер матрицы:
     *
     *     Ne x Nv.
     */
    A->rows = mesh->n_edges;
    A->cols = mesh->n_nodes;

    /*
     * Используем calloc(), потому что подавляющее большинство
     * элементов матрицы должно быть равно нулю.
     *
     * Несмотря на то, что матрица разреженная, пока храним ее
     * в плотном виде для простоты дальнейших операций.
     */
    A->data = calloc(A->rows * A->cols, sizeof(double));

    if (A->data == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate incidence matrix.\n");
        incidence_free(A);
        return -1;
    }

    /*
     * Для каждого глобального ребра известны два узла:
     *
     *     a -> b.
     *
     * Глобальная ориентация уже была задана в topology.c.
     *
     * Заполняем:
     *
     *     A[e, a] = -1
     *     A[e, b] = +1.
     */
    for (size_t e = 0; e < mesh->n_edges; ++e)
    {
        int a = mesh->edges[2 * e + 0];
        int b = mesh->edges[2 * e + 1];

        if (a < 0 || b < 0)
        {
            fprintf(stderr, "ERROR: negative node index in edge %zu.\n", e);
            incidence_free(A);
            return -1;
        }

        if ((size_t)a >= mesh->n_nodes || (size_t)b >= mesh->n_nodes)
        {
            fprintf(stderr, "ERROR: node index out of range in edge %zu.\n", e);
            incidence_free(A);
            return -1;
        }

        /*
         * Row-major:
         *
         *     A[e, node] = data[e * Nv + node].
         */
        A->data[e * A->cols + (size_t)a] = -1.0;
        A->data[e * A->cols + (size_t)b] = +1.0;
    }

    return 0;
}


double incidence_check_row_sum(const IncidenceMatrix *A)
{
    if (A == NULL || A->data == NULL)
    {
        return -1.0;
    }

    double max_error = 0.0;

    /*
     * Теоретически каждая строка содержит:
     *
     *     -1 + 1 = 0.
     *
     * Поэтому:
     *
     *     A * 1 = 0.
     */
    for (size_t row = 0; row < A->rows; ++row)
    {
        double sum = 0.0;

        for (size_t col = 0; col < A->cols; ++col)
        {
            sum += A->data[row * A->cols + col];
        }

        double error = fabs(sum);

        if (error > max_error)
        {
            max_error = error;
        }
    }

    return max_error;
}


void incidence_print_info(const IncidenceMatrix *A)
{
    if (A == NULL)
    {
        return;
    }

    printf("\n");
    printf("========================================\n");
    printf("INCIDENCE MATRIX\n");
    printf("========================================\n");
    printf("Rows Ne : %zu\n", A->rows);
    printf("Cols Nv : %zu\n", A->cols);

    /*
     * В каждой строке должно быть ровно два ненулевых элемента:
     *
     *     -1 и +1.
     */
    size_t nonzero = 0;

    if (A->data != NULL)
    {
        for (size_t i = 0; i < A->rows * A->cols; ++i)
        {
            if (A->data[i] != 0.0)
            {
                ++nonzero;
            }
        }
    }

    printf("Nonzero : %zu\n", nonzero);
    printf("Expected: %zu\n", 2 * A->rows);
    printf("========================================\n");
}

/*
 * ============================================================
 * ПОСТРОЕНИЕ:
 *
 *     K = A P A^T
 * ============================================================
 *
 * Для двух глобальных ребер:
 *
 *     e = (a,b)
 *     f = (c,d)
 *
 * при нашей ориентации:
 *
 *     A[e,a] = -1
 *     A[e,b] = +1
 *
 *     A[f,c] = -1
 *     A[f,d] = +1.
 *
 *
 * Тогда:
 *
 *     K_ef =
 *
 *     (-1)(-1) P_ac
 *     + (-1)(+1) P_ad
 *     + (+1)(-1) P_bc
 *     + (+1)(+1) P_bd
 *
 * то есть:
 *
 *     K_ef =
 *
 *       P_ac
 *     - P_ad
 *     - P_bc
 *     + P_bd.
 *
 *
 * Эквивалентная запись:
 *
 *     K_ef =
 *
 *       P[b,d]
 *     - P[b,c]
 *     - P[a,d]
 *     + P[a,c].
 */
int incidence_build_APAT(
    const Mesh *mesh,
    const PotentialMatrix *P,
    Matrix *K,
    int parallel_threads
)
{
    if (mesh == NULL ||
        P == NULL ||
        K == NULL ||
        P->data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    /*
     * Проверяем согласованность:
     *
     *     P должно быть Nv x Nv.
     */
    if (P->n != mesh->n_nodes)
    {
        fprintf(
            stderr,
            "ERROR: incompatible P matrix size in incidence_build_APAT.\n"
        );

        return -1;
    }

    size_t ne = mesh->n_edges;
    size_t nv = mesh->n_nodes;

    if (ne == 0 || nv == 0)
    {
        return -1;
    }

    /*
     * Результат:
     *
     *     K in R^(Ne x Ne).
     */
    if (matrix_alloc(
        K,
        ne,
        ne
    ) != 0)
    {
        return -1;
    }


    /*
     * ========================================================
     * OPENMP
     * ========================================================
     *
     * Матрица:
     *
     *     K = A P A^T
     *
     * симметрична, если P симметрична.
     *
     * Поэтому достаточно вычислить:
     *
     *     f >= e
     *
     * и затем записать:
     *
     *     K_fe = K_ef.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long e_signed = 0; e_signed < (long long)ne; ++e_signed)
    {
        size_t e = (size_t)e_signed;

        /*
         * Концы первого ребра:
         *
         *     a -> b.
         */
        size_t a =
        (size_t)mesh->edges[2 * e + 0];

        size_t b =
        (size_t)mesh->edges[2 * e + 1];


        for (size_t f = e; f < ne; ++f)
        {
            /*
             * Концы второго ребра:
             *
             *     c -> d.
             */
            size_t c =
            (size_t)mesh->edges[2 * f + 0];

            size_t d =
            (size_t)mesh->edges[2 * f + 1];


            /*
             * Формула:
             *
             *     K_ef =
             *
             *       P[a,c]
             *     - P[a,d]
             *     - P[b,c]
             *     + P[b,d].
             */
            double value =
            P->data[a * nv + c]
            - P->data[a * nv + d]
            - P->data[b * nv + c]
            + P->data[b * nv + d];


            K->data[e * ne + f] = value;
            K->data[f * ne + e] = value;
        }
    }

    return 0;
}

/*
 * ============================================================
 * B = P A^T
 * ============================================================
 *
 * Размер:
 *
 *     B : Nv x Ne.
 *
 *
 * Для глобально ориентированного ребра:
 *
 *     e = (a -> b)
 *
 * матрица инцидентности содержит:
 *
 *     A[e,a] = -1
 *     A[e,b] = +1.
 *
 *
 * Поэтому:
 *
 *     B[j,e]
 *
 *     = sum_k P[j,k] A^T[k,e]
 *
 *     = sum_k P[j,k] A[e,k]
 *
 *     = -P[j,a] + P[j,b].
 *
 *
 * То есть:
 *
 *     B[j,e] = P[j,b] - P[j,a].
 *
 *
 * Это позволяет не выполнять обычное плотное
 * матричное умножение:
 *
 *     P * A^T.
 *
 * У каждого ребра в A только два ненулевых элемента,
 * поэтому вычисление имеет сложность:
 *
 *     O(Nv * Ne).
 */
int incidence_build_PAT(
    const Mesh *mesh,
    const PotentialMatrix *P,
    Matrix *PAT,
    int parallel_threads
)
{
    if (mesh == NULL ||
        P == NULL ||
        PAT == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (mesh->edges == NULL ||
        P->data == NULL)
    {
        fprintf(
            stderr,
            "ERROR: invalid input in incidence_build_PAT().\n"
        );

        return -1;
    }

    if (P->n != mesh->n_nodes)
    {
        fprintf(
            stderr,
            "ERROR: P size does not match mesh nodes in "
            "incidence_build_PAT().\n"
        );

        return -1;
    }


    /*
     * PAT имеет размер:
     *
     *     Nv x Ne.
     */
    if (matrix_alloc(
        PAT,
        mesh->n_nodes,
        mesh->n_edges
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate P A^T matrix.\n"
        );

        return -1;
    }


    size_t nv =
    mesh->n_nodes;

    size_t ne =
    mesh->n_edges;


    /*
     * Параллелим по узлам.
     *
     * Каждая строка PAT[j,:] записывается
     * только одним потоком.
     */
    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long j_signed = 0;
         j_signed < (long long)nv;
    ++j_signed)
         {
             size_t j =
             (size_t)j_signed;

             const double *p_row =
             &P->data[j * nv];

             double *pat_row =
             &PAT->data[j * ne];


             for (size_t e = 0;
                  e < ne;
             ++e)
                  {
                      /*
                       * Глобальная ориентация ребра уже была
                       * определена в topology.c:
                       *
                       *     a -> b.
                       *
                       * mesh->edges хранит:
                       *
                       *     [a0,b0,a1,b1,...].
                       */
                      int a =
                      mesh->edges[2 * e + 0];

                      int b =
                      mesh->edges[2 * e + 1];


                      if (a < 0 ||
                          b < 0 ||
                          (size_t)a >= nv ||
                          (size_t)b >= nv)
                      {
                          /*
                           * Такая ситуация означает поврежденную
                           * топологию.
                           */
                          pat_row[e] =
                          NAN;

                          continue;
                      }


                      /*
                       *     (P A^T)[j,e]
                       *
                       *     = P[j,b] - P[j,a].
                       */
                      pat_row[e] =
                      p_row[(size_t)b]
                      - p_row[(size_t)a];
                  }
         }


         /*
          * Проверяем, что параллельный цикл
          * не обнаружил поврежденные ребра.
          */
         for (size_t j = 0;
              j < nv;
    ++j)
              {
                  for (size_t e = 0;
                       e < ne;
                  ++e)
                       {
                           double value =
                           PAT->data[j * ne + e];

                           if (!isfinite(value))
                           {
                               fprintf(
                                   stderr,
                                   "ERROR: invalid value in P A^T "
                                   "at [%zu,%zu].\n",
                                   j,
                                   e
                               );

                               matrix_free(PAT);

                               return -1;
                           }
                       }
              }

              return 0;
}


/*
 * ============================================================
 * B = A S
 * ============================================================
 */
int incidence_build_AS(
    const IncidenceMatrix *A,
    const PotentialMatrix *P,
    Matrix *AS,
    int parallel_threads
)
{
    if (A == NULL ||
        P == NULL ||
        AS == NULL ||
        A->data == NULL ||
        P->data == NULL ||
        parallel_threads < 1)
    {
        return -1;
    }

    if (A->cols != P->n ||
        A->rows == 0 ||
        A->cols == 0)
    {
        fprintf(
            stderr,
            "ERROR: incompatible dimensions in incidence_build_AS().\n"
        );

        return -1;
    }

    size_t ne =
    A->rows;

    size_t nv =
    A->cols;

    if (matrix_alloc(
        AS,
        ne,
        nv
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot allocate A S matrix.\n"
        );

        return -1;
    }

    int error_code =
    0;

    #pragma omp parallel for schedule(static) num_threads(parallel_threads)
    for (long long e_signed = 0;
         e_signed < (long long)ne;
    ++e_signed)
         {
             size_t e =
             (size_t)e_signed;

             size_t node_minus =
             (size_t)-1;

             size_t node_plus =
             (size_t)-1;

             const double *a_row =
             &A->data[e * nv];

             for (size_t j = 0; j < nv; ++j)
             {
                 if (a_row[j] == -1.0)
                 {
                     node_minus =
                     j;
                 }
                 else if (a_row[j] == 1.0)
                 {
                     node_plus =
                     j;
                 }
                 else if (a_row[j] != 0.0)
                 {
                     #pragma omp critical(incidence_as_error)
                     {
                         error_code =
                         -1;
                     }
                 }
             }

             if (node_minus == (size_t)-1 ||
                 node_plus == (size_t)-1)
             {
                 #pragma omp critical(incidence_as_error)
                 {
                     error_code =
                     -1;
                 }

                 continue;
             }

             for (size_t a = 0; a < nv; ++a)
             {
                 double P_aa =
                 P->data[a * nv + a];

                 if (!isfinite(P_aa) ||
                     P_aa <= 0.0)
                 {
                     #pragma omp critical(incidence_as_error)
                     {
                         error_code =
                         -1;
                     }

                     AS->data[e * nv + a] =
                     NAN;

                     continue;
                 }

                 AS->data[e * nv + a] =
                 (
                     P->data[node_plus * nv + a]
                     -
                     P->data[node_minus * nv + a]
                 )
                 / P_aa;
             }
         }

         if (error_code != 0)
         {
             fprintf(
                 stderr,
                 "ERROR: cannot build A S matrix.\n"
             );

             matrix_free(
                 AS
             );

             return -1;
         }

         return 0;
}
