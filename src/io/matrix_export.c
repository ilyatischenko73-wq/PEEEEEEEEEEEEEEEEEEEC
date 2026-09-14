#include "matrix_export.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>


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
        if (buffer[i] == '/' ||
            buffer[i] == '\0')
        {
            char saved =
                buffer[i];

            buffer[i] =
                '\0';

            if (buffer[0] != '\0' &&
                mkdir(buffer, 0755) != 0 &&
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


static int make_path(
    char *buffer,
    size_t buffer_size,
    const char *directory,
    const char *filename
)
{
    if (buffer == NULL ||
        buffer_size == 0 ||
        directory == NULL ||
        filename == NULL)
    {
        return -1;
    }

    int written =
        snprintf(
            buffer,
            buffer_size,
            "%s/%s",
            directory,
            filename
        );

    if (written < 0 ||
        (size_t)written >= buffer_size)
    {
        return -1;
    }

    return 0;
}


int matrix_export_dense_mtx(
    const char *filename,
    const double *data,
    size_t rows,
    size_t cols
)
{
    if (filename == NULL ||
        data == NULL ||
        rows == 0 ||
        cols == 0)
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
            "ERROR: cannot open matrix file '%s'.\n",
            filename
        );

        return -1;
    }

    /*
     * Matrix Market array format stores values column-major.
     *
     * Internal PEEC storage is row-major, therefore we only
     * change traversal order when writing.
     */
    fprintf(
        file,
        "%%%%MatrixMarket matrix array real general\n"
    );

    fprintf(
        file,
        "%% PEEC matrix export\n"
    );

    fprintf(
        file,
        "%zu %zu\n",
        rows,
        cols
    );

    for (size_t col = 0; col < cols; ++col)
    {
        for (size_t row = 0; row < rows; ++row)
        {
            fprintf(
                file,
                "%.17e\n",
                data[row * cols + col]
            );
        }
    }

    fclose(file);

    return 0;
}


int matrix_export_incidence_mtx(
    const char *filename,
    const Mesh *mesh
)
{
    if (filename == NULL ||
        mesh == NULL ||
        mesh->edges == NULL ||
        mesh->n_edges == 0 ||
        mesh->n_nodes == 0)
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

    fprintf(
        file,
        "%%%%MatrixMarket matrix coordinate real general\n"
    );

    fprintf(
        file,
        "%% Incidence matrix A: Ne x Nv\n"
    );

    fprintf(
        file,
        "%zu %zu %zu\n",
        mesh->n_edges,
        mesh->n_nodes,
        2 * mesh->n_edges
    );

    /*
     * Matrix Market indices are 1-based.
     *
     * Global orientation:
     *
     *     a -> b
     *
     *     A[e,a] = -1
     *     A[e,b] = +1.
     */
    for (size_t e = 0; e < mesh->n_edges; ++e)
    {
        int a =
            mesh->edges[2 * e + 0];

        int b =
            mesh->edges[2 * e + 1];

        fprintf(
            file,
            "%zu %d -1.0\n",
            e + 1,
            a + 1
        );

        fprintf(
            file,
            "%zu %d 1.0\n",
            e + 1,
            b + 1
        );
    }

    fclose(file);

    return 0;
}


int matrix_export_resistance_mtx(
    const char *filename,
    const ResistanceMatrix *R
)
{
    if (filename == NULL ||
        R == NULL ||
        R->diagonal == NULL ||
        R->n == 0)
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

    fprintf(
        file,
        "%%%%MatrixMarket matrix coordinate real general\n"
    );

    fprintf(
        file,
        "%% Diagonal PEEC resistance matrix R\n"
    );

    /*
     * Сохраняем все диагональные элементы, включая нули.
     * Для PEC это позволяет явно восстановить размер матрицы.
     */
    fprintf(
        file,
        "%zu %zu %zu\n",
        R->n,
        R->n,
        R->n
    );

    for (size_t i = 0; i < R->n; ++i)
    {
        fprintf(
            file,
            "%zu %zu %.17e\n",
            i + 1,
            i + 1,
            R->diagonal[i]
        );
    }

    fclose(file);

    return 0;
}


int matrix_export_peec_set(
    const char *directory,
    const Mesh *mesh,
    const InductanceMatrix *L,
    const PotentialMatrix *P,
    const Matrix *PAT,
    const ResistanceMatrix *R
)
{
    if (directory == NULL ||
        mesh == NULL ||
        L == NULL ||
        P == NULL ||
        PAT == NULL ||
        R == NULL)
    {
        return -1;
    }

    if (ensure_directory(
            directory
        ) != 0)
    {
        return -1;
    }

    char path[4096];

    if (make_path(
            path,
            sizeof(path),
            directory,
            "L.mtx"
        ) != 0 ||
        matrix_export_dense_mtx(
            path,
            L->data,
            L->n,
            L->n
        ) != 0)
    {
        return -1;
    }

    if (make_path(
            path,
            sizeof(path),
            directory,
            "A.mtx"
        ) != 0 ||
        matrix_export_incidence_mtx(
            path,
            mesh
        ) != 0)
    {
        return -1;
    }

    if (make_path(
            path,
            sizeof(path),
            directory,
            "P.mtx"
        ) != 0 ||
        matrix_export_dense_mtx(
            path,
            P->data,
            P->n,
            P->n
        ) != 0)
    {
        return -1;
    }

    if (make_path(
            path,
            sizeof(path),
            directory,
            "PAT.mtx"
        ) != 0 ||
        matrix_export_dense_mtx(
            path,
            PAT->data,
            PAT->rows,
            PAT->cols
        ) != 0)
    {
        return -1;
    }

    if (make_path(
            path,
            sizeof(path),
            directory,
            "R.mtx"
        ) != 0 ||
        matrix_export_resistance_mtx(
            path,
            R
        ) != 0)
    {
        return -1;
    }

    printf(
        "PEEC matrices saved : %s\n",
        directory
    );

    return 0;
}
