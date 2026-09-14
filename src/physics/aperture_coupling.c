#include "aperture_coupling.h"

#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * ============================================================
 * INIT
 * ============================================================
 */
void aperture_coupling_init(
    ApertureCoupling *coupling
)
{
    if (coupling == NULL)
    {
        return;
    }

    memset(coupling, 0, sizeof(*coupling));
    coupling->initialized = 1;
}


/*
 * ============================================================
 * FREE
 * ============================================================
 */
void aperture_coupling_free(
    ApertureCoupling *coupling
)
{
    if (coupling == NULL)
    {
        return;
    }

    free(coupling->node_mapping);
    free(coupling->cover_edges);
    free(coupling->cut_edges);

    memset(coupling, 0, sizeof(*coupling));
    coupling->initialized = 1;
}


/*
 * ============================================================
 * FIND MAPPED OPEN NODE
 * ============================================================
 *
 * closed_to_open[closed_node] =
 *
 *     mapped open node
 *
 * или:
 *
 *     SIZE_MAX,
 *
 * если closed node не принадлежит aperture boundary.
 */
static int aperture_build_closed_to_open(
    const ApertureCoupling *coupling,
    size_t **closed_to_open_out
)
{
    if (coupling == NULL ||
        closed_to_open_out == NULL ||
        coupling->n_closed_nodes == 0)
    {
        return -1;
    }

    size_t *closed_to_open = malloc(
        coupling->n_closed_nodes * sizeof(size_t)
    );

    if (closed_to_open == NULL)
    {
        return -1;
    }

    for (size_t j = 0; j < coupling->n_closed_nodes; ++j)
    {
        closed_to_open[j] = SIZE_MAX;
    }

    for (size_t k = 0; k < coupling->n_node_mapping; ++k)
    {
        size_t closed_node = coupling->node_mapping[k].closed_node;
        size_t open_node = coupling->node_mapping[k].open_node;

        if (closed_node >= coupling->n_closed_nodes ||
            open_node >= coupling->n_open_nodes)
        {
            free(closed_to_open);
            return -1;
        }

        /*
         * Один closed boundary node не должен иметь
         * два разных отображения.
         */
        if (closed_to_open[closed_node] != SIZE_MAX)
        {
            fprintf(
                stderr,
                "ERROR: duplicate aperture mapping for closed node %zu.\n",
                closed_node
            );

            free(closed_to_open);
            return -1;
        }

        closed_to_open[closed_node] = open_node;
    }

    *closed_to_open_out = closed_to_open;

    return 0;
}


/*
 * ============================================================
 * GEOMETRY DIAGNOSTICS
 * ============================================================
 *
 * These diagnostics are intentionally based only on the actual
 * Mesh arrays used by the C solver.  Therefore all printed node
 * and edge numbers are PEEC local indices, not Gmsh tags and not
 * Python/PyVista indices.
 */

static double aperture_distance3(
    const double *a,
    const double *b
)
{
    double dx =
    a[0] - b[0];

    double dy =
    a[1] - b[1];

    double dz =
    a[2] - b[2];

    return sqrt(
        dx * dx
        +
        dy * dy
        +
        dz * dz
    );
}


static double aperture_mesh_scale(
    const Mesh *mesh
)
{
    if (mesh == NULL ||
        mesh->xyz == NULL ||
        mesh->n_nodes == 0)
    {
        return 1.0;
    }

    double xmin =
    mesh->xyz[0];

    double xmax =
    mesh->xyz[0];

    double ymin =
    mesh->xyz[1];

    double ymax =
    mesh->xyz[1];

    double zmin =
    mesh->xyz[2];

    double zmax =
    mesh->xyz[2];

    for (size_t j = 1;
         j < mesh->n_nodes;
    ++j)
         {
             double x =
             mesh->xyz[3 * j + 0];

             double y =
             mesh->xyz[3 * j + 1];

             double z =
             mesh->xyz[3 * j + 2];

             if (x < xmin)
             {
                 xmin = x;
             }

             if (x > xmax)
             {
                 xmax = x;
             }

             if (y < ymin)
             {
                 ymin = y;
             }

             if (y > ymax)
             {
                 ymax = y;
             }

             if (z < zmin)
             {
                 zmin = z;
             }

             if (z > zmax)
             {
                 zmax = z;
             }
         }

         double sx =
         xmax - xmin;

         double sy =
         ymax - ymin;

         double sz =
         zmax - zmin;

         double scale =
         fmax(
             sx,
              fmax(
                  sy,
                   sz
              )
         );

         if (!isfinite(scale) ||
             scale <= 0.0)
         {
             scale =
             1.0;
         }

         return scale;
}


static int aperture_coupling_print_geometry_diagnostics(
    const ApertureCoupling *coupling,
    const Mesh *closed_mesh,
    const Mesh *open_mesh
)
{
    if (coupling == NULL ||
        closed_mesh == NULL ||
        open_mesh == NULL ||
        closed_mesh->xyz == NULL ||
        closed_mesh->edges == NULL ||
        open_mesh->xyz == NULL)
    {
        return -1;
    }

    double scale =
    fmax(
        aperture_mesh_scale(
            closed_mesh
        ),
         aperture_mesh_scale(
             open_mesh
         )
    );

    /*
     * Geometric equality tolerance.
     *
     * This is intentionally strict because the open and closed
     * aperture-boundary nodes are expected to represent the same
     * physical points.
     */
    double tolerance =
    1.0e-10
    * scale;

    if (tolerance < 1.0e-12)
    {
        tolerance =
        1.0e-12;
    }

    printf("\n");
    printf("============================================================\n");
    printf("APERTURE GEOMETRY DIAGNOSTICS\n");
    printf("============================================================\n");
    printf("Index convention     : PEEC LOCAL INDICES\n");
    printf("Geometry tolerance   : %.9e m\n", tolerance);
    printf("\n");

    /*
     * --------------------------------------------------------
     * NODE MAPPING
     * --------------------------------------------------------
     */
    printf("NODE MAPPING / closed -> open\n");
    printf("------------------------------------------------------------\n");

    double max_mapping_error =
    0.0;

    size_t mapping_error_count =
    0;

    for (size_t k = 0;
         k < coupling->n_node_mapping;
    ++k)
         {
             size_t closed_node =
             coupling->node_mapping[k].closed_node;

             size_t open_node =
             coupling->node_mapping[k].open_node;

             if (closed_node >= closed_mesh->n_nodes ||
                 open_node >= open_mesh->n_nodes)
             {
                 printf(
                     "[%02zu] INVALID INDEX  closed=%zu  open=%zu\n",
                     k,
                     closed_node,
                     open_node
                 );

                 ++mapping_error_count;

                 continue;
             }

             const double *rc =
             &closed_mesh->xyz[
                 3 * closed_node
             ];

             const double *ro =
             &open_mesh->xyz[
                 3 * open_node
             ];

             double error =
             aperture_distance3(
                 rc,
                 ro
             );

             if (error >
                 max_mapping_error)
             {
                 max_mapping_error =
                 error;
             }

             int ok =
             error <= tolerance;

             if (!ok)
             {
                 ++mapping_error_count;
             }

             printf(
                 "[%02zu] closed=%4zu -> open=%4zu   "
                 "error=% .6e m   %s\n",
                 k,
                 closed_node,
                 open_node,
                 error,
                 ok ? "OK" : "ERROR"
             );

             printf(
                 "     closed xyz = % .12e  % .12e  % .12e\n",
                 rc[0],
                 rc[1],
                 rc[2]
             );

             printf(
                 "     open   xyz = % .12e  % .12e  % .12e\n",
                 ro[0],
                 ro[1],
                 ro[2]
             );
         }

         printf("\n");
         printf(
             "Maximum mapping error : %.9e m\n",
             max_mapping_error
         );

         printf(
             "Mapping errors        : %zu / %zu\n",
             mapping_error_count,
             coupling->n_node_mapping
         );

         /*
          * --------------------------------------------------------
          * CUT EDGES
          * --------------------------------------------------------
          */
         printf("\n");
         printf("CUT EDGES USED IN I_gamma = B_gamma I_cover\n");
         printf("------------------------------------------------------------\n");

         size_t cut_error_count =
         0;

         for (size_t k = 0;
              k < coupling->n_cut_edges;
    ++k)
              {
                  const ApertureCutEdge *cut =
                  &coupling->cut_edges[k];

                  if (cut->closed_edge >=
                      closed_mesh->n_edges)
                  {
                      printf(
                          "[%02zu] INVALID EDGE %zu\n",
                          k,
                          cut->closed_edge
                      );

                      ++cut_error_count;

                      continue;
                  }

                  size_t a =
                  (size_t)closed_mesh->edges[
                      2 * cut->closed_edge + 0
                  ];

                  size_t b =
                  (size_t)closed_mesh->edges[
                      2 * cut->closed_edge + 1
                  ];

                  if (a >= closed_mesh->n_nodes ||
                      b >= closed_mesh->n_nodes ||
                      cut->open_node >= open_mesh->n_nodes)
                  {
                      printf(
                          "[%02zu] INVALID NODE INDEX on edge %zu\n",
                          k,
                          cut->closed_edge
                      );

                      ++cut_error_count;

                      continue;
                  }

                  const double *ra =
                  &closed_mesh->xyz[
                      3 * a
                  ];

                  const double *rb =
                  &closed_mesh->xyz[
                      3 * b
                  ];

                  const double *ro =
                  &open_mesh->xyz[
                      3 * cut->open_node
                  ];

                  /*
                   * By construction:
                   *
                   * sign = -1  -> endpoint a is the aperture boundary node
                   * sign = +1  -> endpoint b is the aperture boundary node
                   */
                  size_t boundary_node =
                  cut->sign < 0.0
                  ? a
                  : b;

                  size_t interior_node =
                  cut->sign < 0.0
                  ? b
                  : a;

                  const double *r_boundary =
                  &closed_mesh->xyz[
                      3 * boundary_node
                  ];

                  const double *r_interior =
                  &closed_mesh->xyz[
                      3 * interior_node
                  ];

                  double boundary_error =
                  aperture_distance3(
                      r_boundary,
                      ro
                  );

                  double dx =
                  rb[0] - ra[0];

                  double dy =
                  rb[1] - ra[1];

                  double dz =
                  rb[2] - ra[2];

                  double length =
                  sqrt(
                      dx * dx
                      +
                      dy * dy
                      +
                      dz * dz
                  );

                  int ok =
                  boundary_error <= tolerance
                  &&
                  isfinite(length)
                  &&
                  length > 0.0
                  &&
                  (
                      cut->sign == -1.0
                      ||
                      cut->sign == +1.0
                  );

                  if (!ok)
                  {
                      ++cut_error_count;
                  }

                  printf(
                      "[%02zu] edge=%4zu  %4zu -> %4zu   "
                      "boundary=%4zu  interior=%4zu  "
                      "open=%4zu  sign=%+g   %s\n",
                      k,
                      cut->closed_edge,
                      a,
                      b,
                      boundary_node,
                      interior_node,
                      cut->open_node,
                      cut->sign,
                      ok ? "OK" : "ERROR"
                  );

                  printf(
                      "     a xyz        = % .12e  % .12e  % .12e\n",
                      ra[0],
                      ra[1],
                      ra[2]
                  );

                  printf(
                      "     b xyz        = % .12e  % .12e  % .12e\n",
                      rb[0],
                      rb[1],
                      rb[2]
                  );

                  printf(
                      "     boundary xyz = % .12e  % .12e  % .12e\n",
                      r_boundary[0],
                      r_boundary[1],
                      r_boundary[2]
                  );

                  printf(
                      "     interior xyz = % .12e  % .12e  % .12e\n",
                      r_interior[0],
                      r_interior[1],
                      r_interior[2]
                  );

                  printf(
                      "     open xyz     = % .12e  % .12e  % .12e\n",
                      ro[0],
                      ro[1],
                      ro[2]
                  );

                  printf(
                      "     boundary map error = %.9e m, length = %.9e m\n",
                      boundary_error,
                      length
                  );
              }

              printf("\n");
              printf(
                  "Cut-edge errors       : %zu / %zu\n",
                  cut_error_count,
                  coupling->n_cut_edges
              );

              printf("\n");

              if (mapping_error_count == 0 &&
                  cut_error_count == 0)
              {
                  printf(
                      "APERTURE GEOMETRY CHECK: PASS\n"
                  );
              }
              else
              {
                  printf(
                      "APERTURE GEOMETRY CHECK: FAIL\n"
                  );
              }

              printf("============================================================\n");

              return
              mapping_error_count == 0
              &&
              cut_error_count == 0
              ? 0
              : -1;
}


/*
 * ============================================================
 * VALIDATE
 * ============================================================
 */
int aperture_coupling_validate(
    const ApertureCoupling *coupling
)
{
    if (coupling == NULL ||
        !coupling->initialized ||
        !coupling->built)
    {
        return -1;
    }

    if (coupling->n_closed_nodes == 0 ||
        coupling->n_closed_edges == 0 ||
        coupling->n_open_nodes == 0 ||
        coupling->node_mapping == NULL ||
        coupling->n_node_mapping == 0 ||
        coupling->cover_edges == NULL ||
        coupling->n_cover_edges == 0)
    {
        return -1;
    }

    if (coupling->n_cut_edges > 0 &&
        coupling->cut_edges == NULL)
    {
        return -1;
    }

    if (coupling->n_perimeter_edges +
        coupling->n_cut_edges +
        coupling->n_internal_edges !=
        coupling->n_cover_edges)
    {
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * BUILD
 * ============================================================
 */
int aperture_coupling_build(
    ApertureCoupling *coupling,
    const Mesh *closed_mesh,
    const Mesh *open_mesh,
    const ApertureNodeMapEntry *node_mapping,
    size_t n_node_mapping,
    const size_t *cover_edges,
    size_t n_cover_edges
)
{
    if (coupling == NULL ||
        closed_mesh == NULL ||
        open_mesh == NULL ||
        closed_mesh->edges == NULL ||
        closed_mesh->n_nodes == 0 ||
        closed_mesh->n_edges == 0 ||
        open_mesh->n_nodes == 0 ||
        node_mapping == NULL ||
        n_node_mapping == 0 ||
        cover_edges == NULL ||
        n_cover_edges == 0)
    {
        return -1;
    }

    aperture_coupling_free(coupling);

    coupling->n_closed_nodes = closed_mesh->n_nodes;
    coupling->n_closed_edges = closed_mesh->n_edges;
    coupling->n_open_nodes = open_mesh->n_nodes;

    coupling->node_mapping = malloc(
        n_node_mapping * sizeof(ApertureNodeMapEntry)
    );

    coupling->cover_edges = malloc(
        n_cover_edges * sizeof(size_t)
    );

    if (coupling->node_mapping == NULL ||
        coupling->cover_edges == NULL)
    {
        aperture_coupling_free(coupling);
        return -1;
    }

    memcpy(
        coupling->node_mapping,
        node_mapping,
        n_node_mapping * sizeof(ApertureNodeMapEntry)
    );

    memcpy(
        coupling->cover_edges,
        cover_edges,
        n_cover_edges * sizeof(size_t)
    );

    coupling->n_node_mapping = n_node_mapping;
    coupling->n_cover_edges = n_cover_edges;

    size_t *closed_to_open = NULL;

    if (aperture_build_closed_to_open(
        coupling,
        &closed_to_open
    ) != 0)
    {
        aperture_coupling_free(coupling);
        return -1;
    }

    /*
     * Сначала только считаем число cut edges.
     */
    size_t n_cut = 0;
    size_t n_perimeter = 0;
    size_t n_internal = 0;

    for (size_t k = 0; k < n_cover_edges; ++k)
    {
        size_t edge = cover_edges[k];

        if (edge >= closed_mesh->n_edges)
        {
            fprintf(
                stderr,
                "ERROR: cover edge %zu is outside closed mesh.\n",
                edge
            );

            free(closed_to_open);
            aperture_coupling_free(coupling);
            return -1;
        }

        size_t a = (size_t)closed_mesh->edges[2 * edge + 0];
        size_t b = (size_t)closed_mesh->edges[2 * edge + 1];

        if (a >= closed_mesh->n_nodes ||
            b >= closed_mesh->n_nodes)
        {
            free(closed_to_open);
            aperture_coupling_free(coupling);
            return -1;
        }

        int a_boundary = closed_to_open[a] != SIZE_MAX;
        int b_boundary = closed_to_open[b] != SIZE_MAX;

        if (a_boundary && b_boundary)
        {
            ++n_perimeter;
        }
        else if (!a_boundary && !b_boundary)
        {
            ++n_internal;
        }
        else
        {
            ++n_cut;
        }
    }

    if (n_cut > 0)
    {
        coupling->cut_edges = calloc(
            n_cut,
            sizeof(ApertureCutEdge)
        );

        if (coupling->cut_edges == NULL)
        {
            free(closed_to_open);
            aperture_coupling_free(coupling);
            return -1;
        }
    }

    /*
     * Теперь строим sparse B_gamma.
     */
    size_t cut_index = 0;

    for (size_t k = 0; k < n_cover_edges; ++k)
    {
        size_t edge = cover_edges[k];

        size_t a = (size_t)closed_mesh->edges[2 * edge + 0];
        size_t b = (size_t)closed_mesh->edges[2 * edge + 1];

        int a_boundary = closed_to_open[a] != SIZE_MAX;
        int b_boundary = closed_to_open[b] != SIZE_MAX;

        if (a_boundary == b_boundary)
        {
            continue;
        }

        ApertureCutEdge *cut = &coupling->cut_edges[cut_index++];

        cut->closed_edge = edge;

        /*
         * Глобальная ориентация:
         *
         *     a -> b.
         *
         * Если boundary endpoint = a:
         *
         *     B_gamma contribution = -I_e.
         *
         * Если boundary endpoint = b:
         *
         *     B_gamma contribution = +I_e.
         */
        if (a_boundary)
        {
            cut->open_node = closed_to_open[a];
            cut->sign = -1.0;
        }
        else
        {
            cut->open_node = closed_to_open[b];
            cut->sign = +1.0;
        }
    }

    free(closed_to_open);

    coupling->n_cut_edges = n_cut;
    coupling->n_perimeter_edges = n_perimeter;
    coupling->n_internal_edges = n_internal;
    coupling->built = 1;

    if (aperture_coupling_validate(coupling) != 0)
    {
        aperture_coupling_free(coupling);
        return -1;
    }

    /* Invalid coordinate mappings must never feed the second-stage source. */
    if (aperture_coupling_print_geometry_diagnostics(
            coupling, closed_mesh, open_mesh) != 0) {
        aperture_coupling_free(coupling);
        return -1;
    }

    return 0;
}


/*
 * ============================================================
 * APPLY
 * ============================================================
 *
 *     I_gamma = B_gamma I_cover.
 *
 *
 * B_gamma не хранится как dense Nv x Ne.
 *
 * Для каждого cut edge существует ровно один ненулевой вклад:
 *
 *     I_gamma[open_node] += sign * I_edge.
 */
int aperture_coupling_apply(
    const ApertureCoupling *coupling,
    const double *closed_edge_current,
    double *open_node_current
)
{
    if (aperture_coupling_validate(coupling) != 0 ||
        closed_edge_current == NULL ||
        open_node_current == NULL)
    {
        return -1;
    }

    memset(
        open_node_current,
        0,
        coupling->n_open_nodes * sizeof(double)
    );

    for (size_t k = 0; k < coupling->n_cut_edges; ++k)
    {
        const ApertureCutEdge *cut = &coupling->cut_edges[k];

        open_node_current[cut->open_node] +=
        cut->sign * closed_edge_current[cut->closed_edge];
    }

    return 0;
}


/*
 * ============================================================
 * APPLY P I_gamma DIRECTLY
 * ============================================================
 *
 * Так как:
 *
 *     I_gamma[n]
 *
 * ненулев только в узлах, связанных с cut edges,
 * можно не собирать dense B_gamma и не выполнять полное
 * умножение P на Nv-вектор.
 */
int aperture_coupling_apply_potential_source(
    const ApertureCoupling *coupling,
    const PotentialMatrix *P_open,
    const double *closed_edge_current,
    double *potential_source
)
{
    if (aperture_coupling_validate(coupling) != 0 ||
        P_open == NULL ||
        P_open->data == NULL ||
        P_open->n != coupling->n_open_nodes ||
        closed_edge_current == NULL ||
        potential_source == NULL)
    {
        return -1;
    }

    size_t nv = coupling->n_open_nodes;

    memset(
        potential_source,
        0,
        nv * sizeof(double)
    );

    for (size_t j = 0; j < nv; ++j)
    {
        double value = 0.0;

        for (size_t k = 0; k < coupling->n_cut_edges; ++k)
        {
            const ApertureCutEdge *cut = &coupling->cut_edges[k];

            double I_edge = closed_edge_current[cut->closed_edge];
            double P_jn = P_open->data[j * nv + cut->open_node];

            value += P_jn * cut->sign * I_edge;
        }

        potential_source[j] = value;
    }

    return 0;
}


/*
 * ============================================================
 * SMALL JSON PARSER
 * ============================================================
 *
 * Нужен только для двух полей:
 *
 *     "node_mapping"
 *     "cover_edges".
 *
 *
 * Это не универсальный JSON parser.
 */
static char *aperture_read_text_file(
    const char *filename
)
{
    if (filename == NULL)
    {
        return NULL;
    }

    FILE *file = fopen(filename, "rb");

    if (file == NULL)
    {
        return NULL;
    }

    if (fseek(file, 0, SEEK_END) != 0)
    {
        fclose(file);
        return NULL;
    }

    long size = ftell(file);

    if (size < 0)
    {
        fclose(file);
        return NULL;
    }

    if (fseek(file, 0, SEEK_SET) != 0)
    {
        fclose(file);
        return NULL;
    }

    char *text = malloc((size_t)size + 1);

    if (text == NULL)
    {
        fclose(file);
        return NULL;
    }

    size_t read_count = fread(text, 1, (size_t)size, file);
    fclose(file);

    if (read_count != (size_t)size)
    {
        free(text);
        return NULL;
    }

    text[size] = '\0';

    return text;
}


static const char *aperture_skip_space(
    const char *p
)
{
    while (p != NULL &&
        *p != '\0' &&
        isspace((unsigned char)*p))
    {
        ++p;
    }

    return p;
}


static const char *aperture_find_array(
    const char *text,
    const char *key
)
{
    if (text == NULL ||
        key == NULL)
    {
        return NULL;
    }

    char pattern[128];

    int written = snprintf(
        pattern,
        sizeof(pattern),
                           "\"%s\"",
                           key
    );

    if (written <= 0 ||
        (size_t)written >= sizeof(pattern))
    {
        return NULL;
    }

    /*
     * В JSON могут существовать одноимённые поля
     * на разных уровнях вложенности.
     *
     * Например:
     *
     *     "diagnostics": {
     *         "n_cover_edges": 24
     *     },
     *
     *     "cover_edges": [
     *         ...
     *     ]
     *
     * или даже старый вариант:
     *
     *     "diagnostics": {
     *         "cover_edges": 24
     *     },
     *
     *     "cover_edges": [
     *         ...
     *     ]
     *
     * Поэтому недостаточно взять первое вхождение ключа.
     * Нужно найти именно то вхождение, после которого
     * значение является JSON-массивом.
     */
    const char *search =
    text;

    while (1)
    {
        const char *p =
        strstr(
            search,
            pattern
        );

        if (p == NULL)
        {
            return NULL;
        }

        const char *colon =
        strchr(
            p + written,
            ':'
        );

        if (colon == NULL)
        {
            return NULL;
        }

        const char *value =
        aperture_skip_space(
            colon + 1
        );

        if (value != NULL &&
            *value == '[')
        {
            return value;
        }

        /*
         * Найдено одноимённое поле, но оно не является
         * массивом. Продолжаем поиск дальше по файлу.
         */
        search =
        p + written;
    }
}


static const char *aperture_find_matching(
    const char *begin,
    char open_char,
    char close_char
)
{
    if (begin == NULL ||
        *begin != open_char)
    {
        return NULL;
    }

    int depth = 0;
    int in_string = 0;
    int escaped = 0;

    for (const char *p = begin; *p != '\0'; ++p)
    {
        if (in_string)
        {
            if (escaped)
            {
                escaped = 0;
                continue;
            }

            if (*p == '\\')
            {
                escaped = 1;
            }
            else if (*p == '"')
            {
                in_string = 0;
            }

            continue;
        }

        if (*p == '"')
        {
            in_string = 1;
            continue;
        }

        if (*p == open_char)
        {
            ++depth;
        }
        else if (*p == close_char)
        {
            --depth;

            if (depth == 0)
            {
                return p;
            }
        }
    }

    return NULL;
}


static int aperture_parse_size_t_value(
    const char *begin,
    const char *end,
    const char *key,
    size_t *value
)
{
    if (begin == NULL ||
        end == NULL ||
        key == NULL ||
        value == NULL ||
        begin >= end)
    {
        return -1;
    }

    char pattern[128];

    int written = snprintf(
        pattern,
        sizeof(pattern),
                           "\"%s\"",
                           key
    );

    if (written <= 0 ||
        (size_t)written >= sizeof(pattern))
    {
        return -1;
    }

    const char *p = strstr(begin, pattern);

    if (p == NULL ||
        p >= end)
    {
        return -1;
    }

    p = strchr(p, ':');

    if (p == NULL ||
        p >= end)
    {
        return -1;
    }

    p = aperture_skip_space(p + 1);

    if (p == NULL ||
        p >= end)
    {
        return -1;
    }

    errno = 0;

    char *number_end = NULL;

    unsigned long long parsed =
    strtoull(p, &number_end, 10);

    if (errno != 0 ||
        number_end == p ||
        number_end > end ||
        parsed > SIZE_MAX)
    {
        return -1;
    }

    *value = (size_t)parsed;

    return 0;
}


static int aperture_parse_node_mapping(
    const char *text,
    ApertureNodeMapEntry **mapping_out,
    size_t *count_out
)
{
    if (text == NULL ||
        mapping_out == NULL ||
        count_out == NULL)
    {
        return -1;
    }

    const char *array_begin =
    aperture_find_array(text, "node_mapping");

    if (array_begin == NULL)
    {
        return -1;
    }

    const char *array_end =
    aperture_find_matching(array_begin, '[', ']');

    if (array_end == NULL)
    {
        return -1;
    }

    size_t capacity = 16;
    size_t count = 0;

    ApertureNodeMapEntry *mapping = malloc(
        capacity * sizeof(ApertureNodeMapEntry)
    );

    if (mapping == NULL)
    {
        return -1;
    }

    const char *p = array_begin + 1;

    while (p < array_end)
    {
        p = aperture_skip_space(p);

        if (p >= array_end)
        {
            break;
        }

        if (*p == ',')
        {
            ++p;
            continue;
        }

        if (*p != '{')
        {
            free(mapping);
            return -1;
        }

        const char *object_end =
        aperture_find_matching(p, '{', '}');

        if (object_end == NULL ||
            object_end > array_end)
        {
            free(mapping);
            return -1;
        }

        size_t closed_node = 0;
        size_t open_node = 0;

        if (aperture_parse_size_t_value(
            p,
            object_end,
            "closed_node",
            &closed_node
        ) != 0 ||
        aperture_parse_size_t_value(
            p,
            object_end,
            "open_node",
            &open_node
        ) != 0)
        {
            free(mapping);
            return -1;
        }

        if (count == capacity)
        {
            size_t new_capacity = 2 * capacity;

            ApertureNodeMapEntry *new_mapping = realloc(
                mapping,
                new_capacity * sizeof(ApertureNodeMapEntry)
            );

            if (new_mapping == NULL)
            {
                free(mapping);
                return -1;
            }

            mapping = new_mapping;
            capacity = new_capacity;
        }

        mapping[count].closed_node = closed_node;
        mapping[count].open_node = open_node;
        ++count;

        p = object_end + 1;
    }

    if (count == 0)
    {
        free(mapping);
        return -1;
    }

    *mapping_out = mapping;
    *count_out = count;

    return 0;
}


static int aperture_parse_cover_edges(
    const char *text,
    size_t **edges_out,
    size_t *count_out
)
{
    if (text == NULL ||
        edges_out == NULL ||
        count_out == NULL)
    {
        return -1;
    }

    const char *array_begin =
    aperture_find_array(text, "cover_edges");

    if (array_begin == NULL)
    {
        return -1;
    }

    const char *array_end =
    aperture_find_matching(array_begin, '[', ']');

    if (array_end == NULL)
    {
        return -1;
    }

    size_t capacity = 32;
    size_t count = 0;

    size_t *edges = malloc(
        capacity * sizeof(size_t)
    );

    if (edges == NULL)
    {
        return -1;
    }

    const char *p = array_begin + 1;

    while (p < array_end)
    {
        p = aperture_skip_space(p);

        if (p >= array_end)
        {
            break;
        }

        if (*p == ',')
        {
            ++p;
            continue;
        }

        errno = 0;

        char *number_end = NULL;

        unsigned long long parsed =
        strtoull(p, &number_end, 10);

        if (errno != 0 ||
            number_end == p ||
            number_end > array_end ||
            parsed > SIZE_MAX)
        {
            free(edges);
            return -1;
        }

        if (count == capacity)
        {
            size_t new_capacity = 2 * capacity;

            size_t *new_edges = realloc(
                edges,
                new_capacity * sizeof(size_t)
            );

            if (new_edges == NULL)
            {
                free(edges);
                return -1;
            }

            edges = new_edges;
            capacity = new_capacity;
        }

        edges[count++] = (size_t)parsed;
        p = number_end;
    }

    if (count == 0)
    {
        free(edges);
        return -1;
    }

    *edges_out = edges;
    *count_out = count;

    return 0;
}


/*
 * ============================================================
 * LOAD JSON
 * ============================================================
 */
int aperture_coupling_load_json(
    ApertureCoupling *coupling,
    const Mesh *closed_mesh,
    const Mesh *open_mesh,
    const char *filename
)
{
    if (coupling == NULL ||
        closed_mesh == NULL ||
        open_mesh == NULL ||
        filename == NULL)
    {
        return -1;
    }

    char *text =
    aperture_read_text_file(filename);

    if (text == NULL)
    {
        fprintf(
            stderr,
            "ERROR: cannot read aperture map file: %s\n",
            filename
        );

        return -1;
    }

    ApertureNodeMapEntry *mapping = NULL;
    size_t n_mapping = 0;

    size_t *cover_edges = NULL;
    size_t n_cover_edges = 0;

    if (aperture_parse_node_mapping(
        text,
        &mapping,
        &n_mapping
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot parse node_mapping in aperture map: %s\n",
            filename
        );

        free(text);
        return -1;
    }

    if (aperture_parse_cover_edges(
        text,
        &cover_edges,
        &n_cover_edges
    ) != 0)
    {
        fprintf(
            stderr,
            "ERROR: cannot parse cover_edges in aperture map: %s\n",
            filename
        );

        free(mapping);
        free(text);

        return -1;
    }

    int result =
    aperture_coupling_build(
        coupling,
        closed_mesh,
        open_mesh,
        mapping,
        n_mapping,
        cover_edges,
        n_cover_edges
    );

    free(cover_edges);
    free(mapping);
    free(text);

    return result;
}


/*
 * ============================================================
 * PRINT
 * ============================================================
 */
void aperture_coupling_print_info(
    const ApertureCoupling *coupling
)
{
    if (coupling == NULL)
    {
        return;
    }

    printf("\n");
    printf("============================================================\n");
    printf("APERTURE COUPLING\n");
    printf("============================================================\n");

    printf(
        "Closed nodes       : %zu\n",
        coupling->n_closed_nodes
    );

    printf(
        "Closed edges       : %zu\n",
        coupling->n_closed_edges
    );

    printf(
        "Open nodes         : %zu\n",
        coupling->n_open_nodes
    );

    printf(
        "Boundary mappings  : %zu\n",
        coupling->n_node_mapping
    );

    printf(
        "Cover edges        : %zu\n",
        coupling->n_cover_edges
    );

    printf(
        "Perimeter edges    : %zu\n",
        coupling->n_perimeter_edges
    );

    printf(
        "Cut edges          : %zu\n",
        coupling->n_cut_edges
    );

    printf(
        "Internal edges     : %zu\n",
        coupling->n_internal_edges
    );

    printf(
        "Transfer           : I_gamma = B_gamma I_cover\n"
    );

    printf(
        "B_gamma storage    : sparse cut-edge form\n"
    );

    printf("============================================================\n");
}
