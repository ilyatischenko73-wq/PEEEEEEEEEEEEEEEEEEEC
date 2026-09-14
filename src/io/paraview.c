#include "paraview.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>


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

    size_t n =
    strlen(directory);

    if (n >= sizeof(buffer))
    {
        return -1;
    }

    memcpy(
        buffer,
        directory,
        n + 1
    );

    for (size_t i = 1; i <= n; ++i)
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
                return -1;
            }

            buffer[i] =
            saved;
        }
    }

    return 0;
}


static int make_filename(
    char *buffer,
    size_t size,
    const char *directory,
    const char *prefix,
    const char *kind,
    long frame
)
{
    if (buffer == NULL ||
        directory == NULL ||
        prefix == NULL ||
        kind == NULL)
    {
        return -1;
    }

    int written;

    if (frame >= 0)
    {
        written =
        snprintf(
            buffer,
            size,
            "%s/%s_%s_%06ld.vtk",
            directory,
            prefix,
            kind,
            frame
        );
    }
    else
    {
        written =
        snprintf(
            buffer,
            size,
            "%s/%s_%s.vtk",
            directory,
            prefix,
            kind
        );
    }

    if (written < 0 ||
        (size_t)written >= size)
    {
        return -1;
    }

    return 0;
}


static int vtk_write_surface_header(
    FILE *file,
    const Mesh *mesh
)
{
    if (file == NULL ||
        mesh == NULL ||
        mesh->xyz == NULL ||
        mesh->quads == NULL)
    {
        return -1;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "PEEC surface state\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");

    fprintf(
        file,
        "POINTS %zu double\n",
        mesh->n_nodes
    );

    for (size_t j = 0; j < mesh->n_nodes; ++j)
    {
        const double *r =
        &mesh->xyz[3 * j];

        fprintf(
            file,
            "%.17e %.17e %.17e\n",
            r[0],
            r[1],
            r[2]
        );
    }

    fprintf(
        file,
        "POLYGONS %zu %zu\n",
        mesh->n_quads,
        5 * mesh->n_quads
    );

    for (size_t q = 0; q < mesh->n_quads; ++q)
    {
        fprintf(
            file,
            "4 %d %d %d %d\n",
            mesh->quads[4 * q + 0],
            mesh->quads[4 * q + 1],
            mesh->quads[4 * q + 2],
            mesh->quads[4 * q + 3]
        );
    }

    return 0;
}


static int vtk_write_edges_header(
    FILE *file,
    const Mesh *mesh
)
{
    if (file == NULL ||
        mesh == NULL ||
        mesh->xyz == NULL ||
        mesh->edges == NULL)
    {
        return -1;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "PEEC edge current state\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");

    fprintf(
        file,
        "POINTS %zu double\n",
        mesh->n_nodes
    );

    for (size_t j = 0; j < mesh->n_nodes; ++j)
    {
        const double *r =
        &mesh->xyz[3 * j];

        fprintf(
            file,
            "%.17e %.17e %.17e\n",
            r[0],
            r[1],
            r[2]
        );
    }

    fprintf(
        file,
        "LINES %zu %zu\n",
        mesh->n_edges,
        3 * mesh->n_edges
    );

    for (size_t e = 0; e < mesh->n_edges; ++e)
    {
        fprintf(
            file,
            "2 %d %d\n",
            mesh->edges[2 * e + 0],
            mesh->edges[2 * e + 1]
        );
    }

    return 0;
}


static double safe_sigma(
    const DualMesh *dual,
    const double *charge,
    size_t j
)
{
    if (dual == NULL ||
        charge == NULL ||
        dual->node_region_areas == NULL ||
        j >= dual->n_nodes ||
        dual->node_region_areas[j] <= 0.0)
    {
        return 0.0;
    }

    return charge[j]
    / dual->node_region_areas[j];
}


static void real_current_density(
    const Mesh *mesh,
    const DualMesh *dual,
    const double *current,
    size_t e,
    double J[3]
)
{
    J[0] = 0.0;
    J[1] = 0.0;
    J[2] = 0.0;

    if (mesh == NULL ||
        dual == NULL ||
        current == NULL ||
        mesh->edge_vectors == NULL ||
        dual->edge_region_areas == NULL ||
        e >= mesh->n_edges ||
        e >= dual->n_edges)
    {
        return;
    }

    double area =
    dual->edge_region_areas[e];

    if (area <= 0.0)
    {
        return;
    }

    /*
     * J = (I/w) e_hat.
     *
     * w = S_e / l
     *
     * l_vec = l e_hat
     *
     * => J_vec = I/S_e * l_vec.
     */
    double factor =
    current[e] / area;

    J[0] =
    factor * mesh->edge_vectors[3 * e + 0];

    J[1] =
    factor * mesh->edge_vectors[3 * e + 1];

    J[2] =
    factor * mesh->edge_vectors[3 * e + 2];
}


static int write_real_surface(
    const char *filename,
    double time,
    const Mesh *mesh,
    const DualMesh *dual,
    const double *V,
    const double *Q
)
{
    FILE *file =
    fopen(
        filename,
        "w"
    );

    if (file == NULL)
    {
        return -1;
    }

    if (vtk_write_surface_header(
        file,
        mesh
    ) != 0)
    {
        fclose(file);
        return -1;
    }

    fprintf(
        file,
        "POINT_DATA %zu\n",
        mesh->n_nodes
    );

    fprintf(file, "SCALARS time double 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (size_t j = 0; j < mesh->n_nodes; ++j)
    {
        fprintf(file, "%.17e\n", time);
    }

    fprintf(file, "SCALARS phi_V double 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (size_t j = 0; j < mesh->n_nodes; ++j)
    {
        fprintf(
            file,
            "%.17e\n",
            V != NULL ? V[j] : 0.0
        );
    }

    if (Q != NULL)
    {
        fprintf(file, "SCALARS charge_C double 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");

        for (size_t j = 0; j < mesh->n_nodes; ++j)
        {
            fprintf(file, "%.17e\n", Q[j]);
        }

        fprintf(file, "SCALARS sigma_C_m2 double 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");

        for (size_t j = 0; j < mesh->n_nodes; ++j)
        {
            fprintf(
                file,
                "%.17e\n",
                safe_sigma(
                    dual,
                    Q,
                    j
                )
            );
        }
    }

    fclose(file);

    return 0;
}


static int write_real_edges(
    const char *filename,
    double time,
    const Mesh *mesh,
    const DualMesh *dual,
    const double *I
)
{
    FILE *file =
    fopen(
        filename,
        "w"
    );

    if (file == NULL)
    {
        return -1;
    }

    if (vtk_write_edges_header(
        file,
        mesh
    ) != 0)
    {
        fclose(file);
        return -1;
    }

    fprintf(
        file,
        "CELL_DATA %zu\n",
        mesh->n_edges
    );

    fprintf(file, "SCALARS time double 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (size_t e = 0; e < mesh->n_edges; ++e)
    {
        fprintf(file, "%.17e\n", time);
    }

    fprintf(file, "SCALARS current_A double 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (size_t e = 0; e < mesh->n_edges; ++e)
    {
        fprintf(
            file,
            "%.17e\n",
            I != NULL ? I[e] : 0.0
        );
    }

    fprintf(file, "VECTORS J_A_m double\n");

    for (size_t e = 0; e < mesh->n_edges; ++e)
    {
        double J[3];

        real_current_density(
            mesh,
            dual,
            I,
            e,
            J
        );

        fprintf(
            file,
            "%.17e %.17e %.17e\n",
            J[0],
            J[1],
            J[2]
        );
    }

    fclose(file);

    return 0;
}


int paraview_write_real_state(
    const char *directory,
    const char *prefix,
    size_t frame,
    double time,
    const Mesh *mesh,
    const DualMesh *dual,
    const double *node_voltage,
    const double *node_charge,
    const double *edge_current
)
{
    if (directory == NULL ||
        prefix == NULL ||
        mesh == NULL ||
        dual == NULL)
    {
        return -1;
    }

    if (ensure_directory(
        directory
    ) != 0)
    {
        return -1;
    }

    char surface_file[4096];
    char edge_file[4096];

    if (make_filename(
        surface_file,
        sizeof(surface_file),
                      directory,
                      prefix,
                      "surface",
                      (long)frame
    ) != 0 ||
    make_filename(
        edge_file,
        sizeof(edge_file),
                  directory,
                  prefix,
                  "edges",
                  (long)frame
    ) != 0)
    {
        return -1;
    }

    if (write_real_surface(
        surface_file,
        time,
        mesh,
        dual,
        node_voltage,
        node_charge
    ) != 0)
    {
        return -1;
    }

    if (write_real_edges(
        edge_file,
        time,
        mesh,
        dual,
        edge_current
    ) != 0)
    {
        return -1;
    }

    return 0;
}


static double complex_abs_value(
    Complex z
)
{
    return hypot(z.re, z.im);
}


static int write_complex_surface(
    const char *filename,
    const Mesh *mesh,
    const DualMesh *dual,
    const Complex *V,
    const Complex *Q
)
{
    FILE *file =
    fopen(
        filename,
        "w"
    );

    if (file == NULL)
    {
        return -1;
    }

    if (vtk_write_surface_header(
        file,
        mesh
    ) != 0)
    {
        fclose(file);
        return -1;
    }

    fprintf(file, "POINT_DATA %zu\n", mesh->n_nodes);

    const char *names[3] = {
        "phi_re_V",
        "phi_im_V",
        "phi_abs_V"
    };

    for (int part = 0; part < 3; ++part)
    {
        fprintf(
            file,
            "SCALARS %s double 1\n",
            names[part]
        );

        fprintf(file, "LOOKUP_TABLE default\n");

        for (size_t j = 0; j < mesh->n_nodes; ++j)
        {
            double value = 0.0;

            if (V != NULL)
            {
                if (part == 0)
                {
                    value = V[j].re;
                }
                else if (part == 1)
                {
                    value = V[j].im;
                }
                else
                {
                    value = complex_abs_value(V[j]);
                }
            }

            fprintf(file, "%.17e\n", value);
        }
    }

    if (Q != NULL)
    {
        const char *q_names[3] = {
            "charge_re_C",
            "charge_im_C",
            "charge_abs_C"
        };

        const char *s_names[3] = {
            "sigma_re_C_m2",
            "sigma_im_C_m2",
            "sigma_abs_C_m2"
        };

        for (int part = 0; part < 3; ++part)
        {
            fprintf(
                file,
                "SCALARS %s double 1\n",
                q_names[part]
            );

            fprintf(file, "LOOKUP_TABLE default\n");

            for (size_t j = 0; j < mesh->n_nodes; ++j)
            {
                double q;

                if (part == 0)
                {
                    q = Q[j].re;
                }
                else if (part == 1)
                {
                    q = Q[j].im;
                }
                else
                {
                    q = complex_abs_value(Q[j]);
                }

                fprintf(file, "%.17e\n", q);
            }

            fprintf(
                file,
                "SCALARS %s double 1\n",
                s_names[part]
            );

            fprintf(file, "LOOKUP_TABLE default\n");

            for (size_t j = 0; j < mesh->n_nodes; ++j)
            {
                double q;

                if (part == 0)
                {
                    q = Q[j].re;
                }
                else if (part == 1)
                {
                    q = Q[j].im;
                }
                else
                {
                    q = complex_abs_value(Q[j]);
                }

                double sigma = 0.0;

                if (dual->node_region_areas != NULL &&
                    j < dual->n_nodes &&
                    dual->node_region_areas[j] > 0.0)
                {
                    sigma =
                    q / dual->node_region_areas[j];
                }

                fprintf(file, "%.17e\n", sigma);
            }
        }
    }

    fclose(file);

    return 0;
}


static int write_complex_edges(
    const char *filename,
    const Mesh *mesh,
    const DualMesh *dual,
    const Complex *I
)
{
    FILE *file =
    fopen(
        filename,
        "w"
    );

    if (file == NULL)
    {
        return -1;
    }

    if (vtk_write_edges_header(
        file,
        mesh
    ) != 0)
    {
        fclose(file);
        return -1;
    }

    fprintf(file, "CELL_DATA %zu\n", mesh->n_edges);

    const char *names[3] = {
        "current_re_A",
        "current_im_A",
        "current_abs_A"
    };

    for (int part = 0; part < 3; ++part)
    {
        fprintf(
            file,
            "SCALARS %s double 1\n",
            names[part]
        );

        fprintf(file, "LOOKUP_TABLE default\n");

        for (size_t e = 0; e < mesh->n_edges; ++e)
        {
            double value = 0.0;

            if (I != NULL)
            {
                if (part == 0)
                {
                    value = I[e].re;
                }
                else if (part == 1)
                {
                    value = I[e].im;
                }
                else
                {
                    value = complex_abs_value(I[e]);
                }
            }

            fprintf(file, "%.17e\n", value);
        }
    }

    /*
     * Сохраняем вектор J отдельно для real и imag.
     */
    for (int part = 0; part < 2; ++part)
    {
        fprintf(
            file,
            "VECTORS %s double\n",
            part == 0
            ? "J_re_A_m"
            : "J_im_A_m"
        );

        for (size_t e = 0; e < mesh->n_edges; ++e)
        {
            double current = 0.0;

            if (I != NULL)
            {
                current =
                part == 0
                ? I[e].re
                : I[e].im;
            }

            double J[3];

            J[0] = 0.0;
            J[1] = 0.0;
            J[2] = 0.0;

            if (mesh->edge_vectors != NULL &&
                dual->edge_region_areas != NULL &&
                e < dual->n_edges &&
                dual->edge_region_areas[e] > 0.0)
            {
                double factor =
                current / dual->edge_region_areas[e];

                J[0] =
                factor * mesh->edge_vectors[3 * e + 0];

                J[1] =
                factor * mesh->edge_vectors[3 * e + 1];

                J[2] =
                factor * mesh->edge_vectors[3 * e + 2];
            }

            fprintf(
                file,
                "%.17e %.17e %.17e\n",
                J[0],
                J[1],
                J[2]
            );
        }
    }

    fclose(file);

    return 0;
}


int paraview_write_complex_state(
    const char *directory,
    const char *prefix,
    const Mesh *mesh,
    const DualMesh *dual,
    const Complex *node_voltage,
    const Complex *node_charge,
    const Complex *edge_current
)
{
    if (directory == NULL ||
        prefix == NULL ||
        mesh == NULL ||
        dual == NULL)
    {
        return -1;
    }

    if (ensure_directory(
        directory
    ) != 0)
    {
        return -1;
    }

    char surface_file[4096];
    char edge_file[4096];

    if (make_filename(
        surface_file,
        sizeof(surface_file),
                      directory,
                      prefix,
                      "surface",
                      -1
    ) != 0 ||
    make_filename(
        edge_file,
        sizeof(edge_file),
                  directory,
                  prefix,
                  "edges",
                  -1
    ) != 0)
    {
        return -1;
    }

    if (write_complex_surface(
        surface_file,
        mesh,
        dual,
        node_voltage,
        node_charge
    ) != 0)
    {
        return -1;
    }

    if (write_complex_edges(
        edge_file,
        mesh,
        dual,
        edge_current
    ) != 0)
    {
        return -1;
    }

    return 0;
}


void paraview_pvd_init(
    ParaViewPvd *pvd
)
{
    if (pvd == NULL)
    {
        return;
    }

    pvd->file = NULL;
    pvd->opened = 0;
}


int paraview_pvd_begin(
    ParaViewPvd *pvd,
    const char *filename
)
{
    if (pvd == NULL ||
        filename == NULL)
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

    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(
        file,
        "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
    );
    fprintf(file, "  <Collection>\n");

    pvd->file = file;
    pvd->opened = 1;

    return 0;
}


int paraview_pvd_add(
    ParaViewPvd *pvd,
    double time,
    const char *vtk_file
)
{
    if (pvd == NULL ||
        !pvd->opened ||
        pvd->file == NULL ||
        vtk_file == NULL)
    {
        return -1;
    }

    FILE *file =
    (FILE *)pvd->file;

    fprintf(
        file,
        "    <DataSet timestep=\"%.17e\" group=\"\" part=\"0\" file=\"%s\"/>\n",
        time,
        vtk_file
    );

    return 0;
}


void paraview_pvd_end(
    ParaViewPvd *pvd
)
{
    if (pvd == NULL ||
        !pvd->opened ||
        pvd->file == NULL)
    {
        return;
    }

    FILE *file =
    (FILE *)pvd->file;

    fprintf(file, "  </Collection>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);

    pvd->file = NULL;
    pvd->opened = 0;
}
