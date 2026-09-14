#include "topology.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


/*
 * Элемент hash-таблицы уникальных ребер.
 *
 * В hash-таблице хранится канонический ключ ребра:
 *
 * node_a < node_b.
 *
 * Реальная ориентация глобального ребра хранится отдельно
 * в mesh->edges и задается первым появлением ребра.
 */
typedef struct
{
    int node_a;
    int node_b;
    int edge_id;
    int used;

} EdgeHashEntry;


/*
 * Hash для пары индексов узлов.
 *
 * Перед вызовом функции всегда выполняется:
 *
 * a < b.
 *
 * Поэтому одно геометрическое ребро всегда имеет
 * один и тот же ключ независимо от локальной ориентации
 * в соседних четырехугольниках.
 */
static uint64_t hash_edge(int a, int b)
{
    uint64_t x = (uint32_t)a;
    uint64_t y = (uint32_t)b;

    /*
     * Объединяем два 32-битных индекса в 64-битное значение.
     */
    uint64_t value = ((x << 32) ^ y);

    /*
     * Дополнительное перемешивание битов.
     *
     * Это уменьшает кластеризацию ключей
     * в hash-таблице с открытой адресацией.
     */
    value ^= (value >> 33);
    value *= UINT64_C(0xff51afd7ed558ccd);

    value ^= (value >> 33);
    value *= UINT64_C(0xc4ceb9fe1a85ec53);

    value ^= (value >> 33);

    return value;
}


/*
 * Возвращает ближайшую степень двойки >= x.
 *
 * Размер hash-таблицы выбирается степенью двойки,
 * чтобы индекс можно было получать через побитовое AND:
 *
 * slot = hash & (table_size - 1)
 */
static size_t next_power_of_two(size_t x)
{
    size_t value = 1;

    while (value < x)
    {
        value <<= 1;
    }

    return value;
}


/*
 * Евклидова норма трехмерного вектора.
 */
static double vector_norm3(double x, double y, double z)
{
    return sqrt(x * x + y * y + z * z);
}


/*
 * Векторное произведение:
 *
 * result = a x b
 */
static void cross3(const double a[3], const double b[3], double result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}


int mesh_build_topology(Mesh *mesh)
{
    if (mesh == NULL || mesh->xyz == NULL || mesh->quads == NULL)
    {
        return -1;
    }

    const size_t nq = mesh->n_quads;

    /*
     * Максимально возможное количество ребер:
     *
     * 4 ребра на каждую ячейку,
     *
     * если бы ни одно ребро не было общим.
     */
    const size_t max_edges = 4 * nq;

    mesh->edges = malloc(2 * max_edges * sizeof(int));
    mesh->cell_edges = malloc(4 * nq * sizeof(int));
    mesh->cell_edge_signs = malloc(4 * nq * sizeof(int));

    if (mesh->edges == NULL ||
        mesh->cell_edges == NULL ||
        mesh->cell_edge_signs == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate topology arrays.\n");
        return -1;
    }

    /*
     * Размер hash-таблицы берем заметно больше ожидаемого
     * числа ребер, чтобы снизить число коллизий.
     */
    size_t table_size = next_power_of_two(8 * nq);

    EdgeHashEntry *table = calloc(table_size, sizeof(EdgeHashEntry));

    if (table == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate edge hash table.\n");
        return -1;
    }

    size_t ne = 0;

    /*
     * Проходим по всем четырехугольным ячейкам.
     */
    for (size_t cell = 0; cell < nq; ++cell)
    {
        const int *q = &mesh->quads[4 * cell];

        /*
         * Локальная ориентация ребер quadrangle:
         *
         * q0 -> q1
         * q1 -> q2
         * q2 -> q3
         * q3 -> q0
         */
        int local_edges[4][2] = {
            {q[0], q[1]},
            {q[1], q[2]},
            {q[2], q[3]},
            {q[3], q[0]}
        };

        for (int local = 0; local < 4; ++local)
        {
            int local_a = local_edges[local][0];
            int local_b = local_edges[local][1];

            /*
             * Python-совместимая нумерация/ориентация.
             *
             * В Python canonical pair (min,max) используется
             * ТОЛЬКО как ключ поиска одного и того же
             * геометрического ребра.
             *
             * Ориентация самого глобального ребра задается
             * первым встретившимся локальным направлением
             * (local_a -> local_b).
             */
            int key_a = local_a;
            int key_b = local_b;

            if (key_a > key_b)
            {
                int tmp = key_a;
                key_a = key_b;
                key_b = tmp;
            }

            uint64_t hash = hash_edge(key_a, key_b);

            /*
             * table_size является степенью двойки,
             * поэтому вместо % используем побитовое AND.
             */
            size_t slot = ((size_t)hash & (table_size - 1));

            /*
             * Открытая адресация с линейным пробированием.
             */
            for (;;)
            {
                EdgeHashEntry *entry = &table[slot];

                /*
                 * Пустая ячейка означает, что такого ребра
                 * раньше еще не было.
                 */
                if (!entry->used)
                {
                    entry->used = 1;

                    /*
                     * В hash-таблице хранится canonical pair,
                     * чтобы поиск не зависел от ориентации.
                     */
                    entry->node_a = key_a;
                    entry->node_b = key_b;
                    entry->edge_id = (int)ne;

                    /*
                     * Но само глобальное ребро сохраняем ровно
                     * в ориентации его первого появления.
                     * Это повторяет mesh_parser.py.
                     */
                    mesh->edges[2 * ne + 0] = local_a;
                    mesh->edges[2 * ne + 1] = local_b;

                    mesh->cell_edges[4 * cell + local] = (int)ne;
                    mesh->cell_edge_signs[4 * cell + local] = +1;

                    ++ne;

                    break;
                }

                /*
                 * Если геометрическое ребро уже существует,
                 * сравниваем текущее локальное направление
                 * с сохраненной глобальной ориентацией.
                 */
                if (entry->node_a == key_a &&
                    entry->node_b == key_b)
                {
                    int edge_id = entry->edge_id;

                    int stored_a =
                    mesh->edges[2 * edge_id + 0];

                    int stored_b =
                    mesh->edges[2 * edge_id + 1];

                    int sign = 0;

                    if (local_a == stored_a &&
                        local_b == stored_b)
                    {
                        sign = +1;
                    }
                    else if (local_a == stored_b &&
                        local_b == stored_a)
                    {
                        sign = -1;
                    }
                    else
                    {
                        fprintf(
                            stderr,
                            "ERROR: inconsistent edge orientation.\n"
                        );

                        free(table);
                        return -1;
                    }

                    mesh->cell_edges[4 * cell + local] = edge_id;
                    mesh->cell_edge_signs[4 * cell + local] = sign;

                    break;
                }

                /*
                 * Коллизия hash-таблицы:
                 * переходим к следующему слоту.
                 */
                slot = ((slot + 1) & (table_size - 1));
            }
        }
    }

    free(table);

    mesh->n_edges = ne;

    /*
     * Уменьшаем массив ребер до фактического размера.
     */
    int *new_edges = realloc(mesh->edges, 2 * ne * sizeof(int));

    if (new_edges != NULL)
    {
        mesh->edges = new_edges;
    }

    /*
     * Выделяем память для геометрических характеристик ребер.
     */
    mesh->edge_vectors = malloc(3 * ne * sizeof(double));
    mesh->edge_lengths = malloc(ne * sizeof(double));
    mesh->edge_directions = malloc(3 * ne * sizeof(double));
    mesh->edge_centers = malloc(3 * ne * sizeof(double));

    /*
     * Выделяем память для геометрии четырехугольных ячеек.
     */
    mesh->cell_centers = malloc(3 * nq * sizeof(double));
    mesh->cell_areas = malloc(nq * sizeof(double));
    mesh->cell_normals = malloc(3 * nq * sizeof(double));

    if (mesh->edge_vectors == NULL ||
        mesh->edge_lengths == NULL ||
        mesh->edge_directions == NULL ||
        mesh->edge_centers == NULL ||
        mesh->cell_centers == NULL ||
        mesh->cell_areas == NULL ||
        mesh->cell_normals == NULL)
    {
        fprintf(stderr, "ERROR: cannot allocate mesh geometry arrays.\n");
        return -1;
    }

    /*
     * =======================================================
     * ГЕОМЕТРИЯ РЕБЕР
     * =======================================================
     */
    for (size_t e = 0; e < ne; ++e)
    {
        int a = mesh->edges[2 * e + 0];
        int b = mesh->edges[2 * e + 1];

        double ax = mesh->xyz[3 * a + 0];
        double ay = mesh->xyz[3 * a + 1];
        double az = mesh->xyz[3 * a + 2];

        double bx = mesh->xyz[3 * b + 0];
        double by = mesh->xyz[3 * b + 1];
        double bz = mesh->xyz[3 * b + 2];

        /*
         * Вектор глобально ориентированного ребра:
         *
         * d = r_b - r_a.
         */
        double dx = bx - ax;
        double dy = by - ay;
        double dz = bz - az;

        /*
         * Длина ребра:
         *
         * l = |d|.
         */
        double length = vector_norm3(dx, dy, dz);

        if (length <= 0.0)
        {
            fprintf(stderr, "ERROR: zero-length edge %zu.\n", e);
            return -1;
        }

        mesh->edge_vectors[3 * e + 0] = dx;
        mesh->edge_vectors[3 * e + 1] = dy;
        mesh->edge_vectors[3 * e + 2] = dz;

        mesh->edge_lengths[e] = length;

        /*
         * Единичный вектор направления:
         *
         * e_i = d / |d|.
         */
        mesh->edge_directions[3 * e + 0] = dx / length;
        mesh->edge_directions[3 * e + 1] = dy / length;
        mesh->edge_directions[3 * e + 2] = dz / length;

        /*
         * Центр ребра:
         *
         * r_c = (r_a + r_b) / 2.
         */
        mesh->edge_centers[3 * e + 0] = 0.5 * (ax + bx);
        mesh->edge_centers[3 * e + 1] = 0.5 * (ay + by);
        mesh->edge_centers[3 * e + 2] = 0.5 * (az + bz);
    }

    /*
     * =======================================================
     * ГЕОМЕТРИЯ ЧЕТЫРЕХУГОЛЬНЫХ ЯЧЕЕК
     * =======================================================
     *
     * Площадь quadrangle вычисляем как сумму площадей
     * двух треугольников:
     *
     * (q0, q1, q2)
     *
     * и
     *
     * (q0, q2, q3).
     *
     * Это совпадает с принятой ранее Python-реализацией.
     */
    for (size_t cell = 0; cell < nq; ++cell)
    {
        const int *q = &mesh->quads[4 * cell];

        double p[4][3];

        for (int k = 0; k < 4; ++k)
        {
            p[k][0] = mesh->xyz[3 * q[k] + 0];
            p[k][1] = mesh->xyz[3 * q[k] + 1];
            p[k][2] = mesh->xyz[3 * q[k] + 2];
        }

        /*
         * Геометрический центр quadrangle:
         *
         * среднее арифметическое четырех вершин.
         */
        mesh->cell_centers[3 * cell + 0] =
        0.25 * (p[0][0] + p[1][0] + p[2][0] + p[3][0]);

        mesh->cell_centers[3 * cell + 1] =
        0.25 * (p[0][1] + p[1][1] + p[2][1] + p[3][1]);

        mesh->cell_centers[3 * cell + 2] =
        0.25 * (p[0][2] + p[1][2] + p[2][2] + p[3][2]);

        /*
         * Векторы первого треугольника:
         *
         * q0 -> q1
         * q0 -> q2
         */
        double v01[3] = {
            p[1][0] - p[0][0],
            p[1][1] - p[0][1],
            p[1][2] - p[0][2]
        };

        double v02[3] = {
            p[2][0] - p[0][0],
            p[2][1] - p[0][1],
            p[2][2] - p[0][2]
        };

        /*
         * Вектор q0 -> q3 используется для второго треугольника.
         */
        double v03[3] = {
            p[3][0] - p[0][0],
            p[3][1] - p[0][1],
            p[3][2] - p[0][2]
        };

        double cross1[3];
        double cross2[3];

        /*
         * Векторные произведения пропорциональны
         * ориентированным площадям треугольников.
         */
        cross3(v01, v02, cross1);
        cross3(v02, v03, cross2);

        double norm1 = vector_norm3(cross1[0], cross1[1], cross1[2]);
        double norm2 = vector_norm3(cross2[0], cross2[1], cross2[2]);

        /*
         * Площадь каждого треугольника равна половине
         * нормы соответствующего векторного произведения.
         */
        double area = 0.5 * (norm1 + norm2);

        mesh->cell_areas[cell] = area;

        /*
         * Нормаль quadrangle определяем через сумму
         * ориентированных нормалей двух треугольников.
         */
        double nx = cross1[0] + cross2[0];
        double ny = cross1[1] + cross2[1];
        double nz = cross1[2] + cross2[2];

        double normal_length = vector_norm3(nx, ny, nz);

        /*
         * Нормализуем нормаль.
         */
        if (normal_length > 0.0)
        {
            nx /= normal_length;
            ny /= normal_length;
            nz /= normal_length;
        }

        mesh->cell_normals[3 * cell + 0] = nx;
        mesh->cell_normals[3 * cell + 1] = ny;
        mesh->cell_normals[3 * cell + 2] = nz;
    }

    return 0;
}
