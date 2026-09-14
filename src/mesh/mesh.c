#include "mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void mesh_init(Mesh *mesh)
{
    if (mesh == NULL)
    {
        return;
    }

    /*
     * Обнуляем все поля структуры.
     * После этого все указатели равны NULL,
     * а размеры равны нулю.
     */
    memset(mesh, 0, sizeof(*mesh));
}


void mesh_free(Mesh *mesh)
{
    if (mesh == NULL)
    {
        return;
    }

    /*
     * Освобождаем данные исходной сетки.
     */
    free(mesh->xyz);
    free(mesh->quads);

    /*
     * Освобождаем топологические данные.
     */
    free(mesh->edges);
    free(mesh->cell_edges);
    free(mesh->cell_edge_signs);

    /*
     * Освобождаем геометрию ребер.
     */
    free(mesh->edge_vectors);
    free(mesh->edge_lengths);
    free(mesh->edge_directions);
    free(mesh->edge_centers);

    /*
     * Освобождаем геометрию ячеек.
     */
    free(mesh->cell_centers);
    free(mesh->cell_areas);
    free(mesh->cell_normals);

    /*
     * После освобождения снова переводим структуру
     * в корректное пустое состояние.
     */
    mesh_init(mesh);
}


void mesh_print_info(const Mesh *mesh)
{
    if (mesh == NULL)
    {
        return;
    }

    printf("\n");
    printf("========================================\n");
    printf("MESH\n");
    printf("========================================\n");
    printf("Nodes : %zu\n", mesh->n_nodes);
    printf("Edges : %zu\n", mesh->n_edges);
    printf("Quads : %zu\n", mesh->n_quads);
    printf("========================================\n");
}
