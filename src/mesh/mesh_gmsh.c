#include "mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define LINE_BUFFER_SIZE 4096


/*
 * Таблица соответствия:
 *
 * Gmsh node tag -> локальный индекс узла в Mesh.
 *
 * Например:
 *
 * Gmsh tag = 125
 *
 * может соответствовать:
 *
 * Mesh index = 37
 *
 * Внутри вычислительного ядра всегда используем компактные
 * индексы 0 ... n_nodes-1.
 */
typedef struct
{
    size_t size;
    int *data;

} NodeTagMap;


/*
 * Освобождение карты тегов узлов.
 */
static void node_tag_map_free(NodeTagMap *map)
{
    if (map == NULL)
    {
        return;
    }

    free(map->data);

    map->data = NULL;
    map->size = 0;
}


/*
 * Создает плотную таблицу Gmsh tag -> local index.
 *
 * Для наших текущих сеток этого достаточно.
 * Позже при необходимости можно заменить на hash map.
 */
static int node_tag_map_create(NodeTagMap *map, size_t max_tag)
{
    map->size = max_tag + 1;

    map->data = malloc(map->size * sizeof(int));

    if (map->data == NULL)
    {
        return -1;
    }

    /*
     * Значение -1 означает, что данный tag еще не зарегистрирован.
     */
    for (size_t i = 0; i < map->size; ++i)
    {
        map->data[i] = -1;
    }

    return 0;
}


/*
 * Читает секцию:
 *
 * $MeshFormat
 * 4.1 0 8
 * $EndMeshFormat
 *
 * file_type = 0 означает ASCII.
 */
static int parse_mesh_format(FILE *file)
{
    char line[LINE_BUFFER_SIZE];

    if (fgets(line, sizeof(line), file) == NULL)
    {
        return -1;
    }

    double version = 0.0;
    int file_type = 0;
    int data_size = 0;

    if (sscanf(line, "%lf %d %d", &version, &file_type, &data_size) != 3)
    {
        fprintf(stderr, "ERROR: invalid $MeshFormat line.\n");
        return -1;
    }

    /*
     * На первом этапе поддерживаем только ASCII Gmsh.
     */
    if (file_type != 0)
    {
        fprintf(stderr, "ERROR: binary Gmsh is not supported yet.\n");
        return -1;
    }

    /*
     * Сейчас ожидаем формат Gmsh 4.x.
     */
    if (version < 4.0 || version >= 5.0)
    {
        fprintf(stderr, "ERROR: expected Gmsh 4.x, got %.3f.\n", version);
        return -1;
    }

    printf("Gmsh version : %.3f\n", version);

    /*
     * Читаем строку $EndMeshFormat.
     */
    if (fgets(line, sizeof(line), file) == NULL)
    {
        return -1;
    }

    return 0;
}


/*
 * Читает секцию $Nodes формата Gmsh 4.x.
 *
 * В Gmsh узлы сгруппированы по entity blocks.
 * Сначала для блока идут node tags, затем координаты.
 */
static int read_nodes(FILE *file, Mesh *mesh, NodeTagMap *tag_map)
{
    char line[LINE_BUFFER_SIZE];

    size_t n_blocks = 0;
    size_t n_nodes = 0;
    size_t min_tag = 0;
    size_t max_tag = 0;

    /*
     * Заголовок секции $Nodes:
     *
     * numEntityBlocks numNodes minNodeTag maxNodeTag
     */
    if (fgets(line, sizeof(line), file) == NULL)
    {
        return -1;
    }

    if (sscanf(line, "%zu %zu %zu %zu",
        &n_blocks,
        &n_nodes,
        &min_tag,
        &max_tag) != 4)
    {
        fprintf(stderr, "ERROR: invalid $Nodes header.\n");
        return -1;
    }

    /*
     * min_tag пока не используется, но читается для соответствия формату.
     */
    (void)min_tag;

    mesh->n_nodes = n_nodes;

    /*
     * Выделяем память под координаты всех узлов.
     *
     * На один узел приходится три double:
     *
     * x, y, z.
     */
    mesh->xyz = malloc(3 * n_nodes * sizeof(double));

    if (mesh->xyz == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate mesh node coordinates.\n");
        return -1;
    }

    /*
     * Создаем таблицу преобразования Gmsh tag -> local index.
     */
    if (node_tag_map_create(tag_map, max_tag) != 0)
    {
        fprintf(stderr, "ERROR: cannot allocate node tag map.\n");
        return -1;
    }

    size_t local_node_index = 0;

    /*
     * Обрабатываем все entity blocks.
     */
    for (size_t block = 0; block < n_blocks; ++block)
    {
        int entity_dim = 0;
        int entity_tag = 0;
        int parametric = 0;
        size_t n_block_nodes = 0;

        /*
         * Заголовок блока:
         *
         * entityDim entityTag parametric numNodesInBlock
         */
        if (fgets(line, sizeof(line), file) == NULL)
        {
            return -1;
        }

        if (sscanf(line, "%d %d %d %zu",
            &entity_dim,
            &entity_tag,
            &parametric,
            &n_block_nodes) != 4)
        {
            fprintf(stderr, "ERROR: invalid node block header.\n");
            return -1;
        }

        /*
         * На текущем этапе размерность и tag entity
         * непосредственно не нужны.
         */
        (void)entity_dim;
        (void)entity_tag;

        /*
         * Параметрические координаты пока не поддерживаем.
         */
        if (parametric != 0)
        {
            fprintf(stderr, "ERROR: parametric Gmsh nodes are not supported yet.\n");
            return -1;
        }

        /*
         * Во временном массиве сохраняем Gmsh tags узлов данного блока.
         */
        size_t *tags = malloc(n_block_nodes * sizeof(size_t));

        if (tags == NULL)
        {
            fprintf(stderr, "ERROR: cannot allocate temporary node tags.\n");
            return -1;
        }

        /*
         * Сначала Gmsh записывает все node tags блока.
         */
        for (size_t i = 0; i < n_block_nodes; ++i)
        {
            if (fgets(line, sizeof(line), file) == NULL)
            {
                free(tags);
                return -1;
            }

            if (sscanf(line, "%zu", &tags[i]) != 1)
            {
                fprintf(stderr, "ERROR: invalid Gmsh node tag.\n");
                free(tags);
                return -1;
            }
        }

        /*
         * После tags идут координаты в том же порядке.
         */
        for (size_t i = 0; i < n_block_nodes; ++i)
        {
            if (fgets(line, sizeof(line), file) == NULL)
            {
                free(tags);
                return -1;
            }

            double x = 0.0;
            double y = 0.0;
            double z = 0.0;

            if (sscanf(line, "%lf %lf %lf", &x, &y, &z) != 3)
            {
                fprintf(stderr, "ERROR: invalid Gmsh node coordinates.\n");
                free(tags);
                return -1;
            }

            if (local_node_index >= n_nodes)
            {
                fprintf(stderr, "ERROR: number of nodes exceeds header value.\n");
                free(tags);
                return -1;
            }

            /*
             * Сохраняем координаты в компактном массиве Mesh.
             */
            mesh->xyz[3 * local_node_index + 0] = x;
            mesh->xyz[3 * local_node_index + 1] = y;
            mesh->xyz[3 * local_node_index + 2] = z;

            size_t gmsh_tag = tags[i];

            if (gmsh_tag >= tag_map->size)
            {
                fprintf(stderr, "ERROR: Gmsh node tag out of range: %zu.\n", gmsh_tag);
                free(tags);
                return -1;
            }

            /*
             * Запоминаем соответствие:
             *
             * Gmsh tag -> local zero-based index.
             */
            tag_map->data[gmsh_tag] = (int)local_node_index;

            ++local_node_index;
        }

        free(tags);
    }

    /*
     * Проверяем, что количество реально прочитанных узлов
     * совпадает с количеством из заголовка.
     */
    if (local_node_index != n_nodes)
    {
        fprintf(stderr,
                "ERROR: expected %zu nodes, but read %zu.\n",
                n_nodes,
                local_node_index);

        return -1;
    }

    /*
     * Следующая строка должна быть $EndNodes.
     */
    if (fgets(line, sizeof(line), file) == NULL)
    {
        return -1;
    }

    return 0;
}


/*
 * Читает секцию $Elements.
 *
 * В текущей PEEC-модели используем только четырехугольные
 * поверхностные элементы Gmsh type = 3.
 *
 * Остальные элементы пока просто пропускаются.
 */
static int read_elements(FILE *file, Mesh *mesh, const NodeTagMap *tag_map)
{
    char line[LINE_BUFFER_SIZE];

    size_t n_blocks = 0;
    size_t n_elements = 0;
    size_t min_tag = 0;
    size_t max_tag = 0;

    /*
     * Заголовок:
     *
     * numEntityBlocks numElements minElementTag maxElementTag
     */
    if (fgets(line, sizeof(line), file) == NULL)
    {
        return -1;
    }

    if (sscanf(line, "%zu %zu %zu %zu",
        &n_blocks,
        &n_elements,
        &min_tag,
        &max_tag) != 4)
    {
        fprintf(stderr, "ERROR: invalid $Elements header.\n");
        return -1;
    }

    /*
     * Пока эти значения нужны только для чтения заголовка.
     */
    (void)n_elements;
    (void)min_tag;
    (void)max_tag;

    /*
     * Начальная емкость массива четырехугольников.
     *
     * Если ее окажется недостаточно, массив будет увеличиваться realloc().
     */
    size_t quad_capacity = 1024;

    mesh->quads = malloc(4 * quad_capacity * sizeof(int));

    if (mesh->quads == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate quad array.\n");
        return -1;
    }

    mesh->n_quads = 0;

    /*
     * Обрабатываем все element blocks.
     */
    for (size_t block = 0; block < n_blocks; ++block)
    {
        int entity_dim = 0;
        int entity_tag = 0;
        int element_type = 0;
        size_t n_block_elements = 0;

        if (fgets(line, sizeof(line), file) == NULL)
        {
            return -1;
        }

        if (sscanf(line, "%d %d %d %zu",
            &entity_dim,
            &entity_tag,
            &element_type,
            &n_block_elements) != 4)
        {
            fprintf(stderr, "ERROR: invalid Gmsh element block header.\n");
            return -1;
        }

        (void)entity_dim;
        (void)entity_tag;

        /*
         * Gmsh element type 3:
         *
         * четырехузловой quadrangle.
         */
        if (element_type == 3)
        {
            for (size_t i = 0; i < n_block_elements; ++i)
            {
                if (fgets(line, sizeof(line), file) == NULL)
                {
                    return -1;
                }

                size_t element_tag = 0;
                size_t n0_tag = 0;
                size_t n1_tag = 0;
                size_t n2_tag = 0;
                size_t n3_tag = 0;

                /*
                 * Для quad строка имеет вид:
                 *
                 * elementTag node1 node2 node3 node4
                 */
                if (sscanf(line, "%zu %zu %zu %zu %zu",
                    &element_tag,
                    &n0_tag,
                    &n1_tag,
                    &n2_tag,
                    &n3_tag) != 5)
                {
                    fprintf(stderr, "ERROR: invalid Gmsh quadrangle.\n");
                    return -1;
                }

                (void)element_tag;

                /*
                 * Если массив заполнен, удваиваем его емкость.
                 */
                if (mesh->n_quads >= quad_capacity)
                {
                    quad_capacity *= 2;

                    int *new_quads = realloc(
                        mesh->quads,
                        4 * quad_capacity * sizeof(int)
                    );

                    if (new_quads == NULL)
                    {
                        fprintf(stderr, "ERROR: cannot grow quad array.\n");
                        return -1;
                    }

                    mesh->quads = new_quads;
                }

                size_t tags[4] = {
                    n0_tag,
                    n1_tag,
                    n2_tag,
                    n3_tag
                };

                /*
                 * Переводим Gmsh node tags в локальные индексы Mesh.
                 */
                for (int k = 0; k < 4; ++k)
                {
                    if (tags[k] >= tag_map->size)
                    {
                        fprintf(stderr, "ERROR: invalid node tag %zu.\n", tags[k]);
                        return -1;
                    }

                    int local_index = tag_map->data[tags[k]];

                    if (local_index < 0)
                    {
                        fprintf(stderr, "ERROR: unknown node tag %zu.\n", tags[k]);
                        return -1;
                    }

                    mesh->quads[4 * mesh->n_quads + k] = local_index;
                }

                ++mesh->n_quads;
            }
        }
        else
        {
            /*
             * Остальные типы элементов пока не используются.
             *
             * Их строки необходимо прочитать и выбросить,
             * чтобы корректно перейти к следующему блоку.
             */
            for (size_t i = 0; i < n_block_elements; ++i)
            {
                if (fgets(line, sizeof(line), file) == NULL)
                {
                    return -1;
                }
            }
        }
    }

    /*
     * Уменьшаем массив до реального количества quadrangles.
     */
    if (mesh->n_quads > 0)
    {
        int *new_quads = realloc(
            mesh->quads,
            4 * mesh->n_quads * sizeof(int)
        );

        if (new_quads != NULL)
        {
            mesh->quads = new_quads;
        }
    }

    /*
     * Читаем $EndElements.
     */
    if (fgets(line, sizeof(line), file) == NULL)
    {
        return -1;
    }

    return 0;
}


int mesh_load_gmsh(const char *filename, Mesh *mesh)
{
    if (filename == NULL || mesh == NULL)
    {
        return -1;
    }

    FILE *file = fopen(filename, "r");

    if (file == NULL)
    {
        fprintf(stderr, "ERROR: cannot open mesh file: %s\n", filename);
        return -1;
    }

    char line[LINE_BUFFER_SIZE];

    NodeTagMap tag_map = {0, NULL};

    int have_nodes = 0;
    int have_elements = 0;

    /*
     * Последовательно просматриваем секции файла.
     *
     * Неизвестные секции просто пропускаются.
     */
    while (fgets(line, sizeof(line), file) != NULL)
    {
        if (strncmp(line, "$MeshFormat", 11) == 0)
        {
            if (parse_mesh_format(file) != 0)
            {
                goto error;
            }
        }
        else if (strncmp(line, "$Nodes", 6) == 0)
        {
            if (read_nodes(file, mesh, &tag_map) != 0)
            {
                goto error;
            }

            have_nodes = 1;
        }
        else if (strncmp(line, "$Elements", 9) == 0)
        {
            /*
             * Элементы ссылаются на node tags,
             * поэтому узлы должны быть прочитаны раньше.
             */
            if (!have_nodes)
            {
                fprintf(stderr, "ERROR: $Elements appeared before $Nodes.\n");
                goto error;
            }

            if (read_elements(file, mesh, &tag_map) != 0)
            {
                goto error;
            }

            have_elements = 1;
        }
    }

    fclose(file);
    node_tag_map_free(&tag_map);

    if (!have_nodes || !have_elements)
    {
        fprintf(stderr, "ERROR: mesh has no usable nodes or elements.\n");
        return -1;
    }

    if (mesh->n_quads == 0)
    {
        fprintf(stderr, "ERROR: mesh has no 4-node quadrangles.\n");
        return -1;
    }

    return 0;


    error:

    fclose(file);
    node_tag_map_free(&tag_map);

    return -1;
}
