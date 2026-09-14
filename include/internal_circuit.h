#ifndef PEEC_INTERNAL_CIRCUIT_H
#define PEEC_INTERNAL_CIRCUIT_H

#include "mesh.h"

#include <stddef.h>


/*
 * ============================================================
 * ТИП ВНУТРЕННЕЙ ЦЕПИ
 * ============================================================
 *
 * Пока реализован только один элемент:
 *
 *     резистивный шунт между двумя узлами корпуса.
 *
 *
 * В дальнейшем сюда можно добавить:
 *
 *     INTERNAL_CIRCUIT_RL
 *     INTERNAL_CIRCUIT_RLC
 *     INTERNAL_CIRCUIT_IMPEDANCE
 *
 * без изменения интерфейса transient solver.
 */
typedef enum
{
    INTERNAL_CIRCUIT_RESISTIVE_SHUNT = 0

} InternalCircuitType;


/*
 * ============================================================
 * ВНУТРЕННИЙ ШУНТ
 * ============================================================
 *
 * Физически:
 *
 *                I_shunt
 *             A --------> B
 *                   R
 *
 *
 * Ориентация ветви:
 *
 *     A -> B.
 *
 *
 * Вектор инцидентности шунта:
 *
 *     b_R[A] = -1
 *     b_R[B] = +1.
 *
 *
 * Поэтому:
 *
 *     b_R^T V = V_B - V_A.
 *
 *
 * Уравнение ветви:
 *
 *     b_R^T V + R I_shunt = 0.
 *
 *
 * То есть:
 *
 *     V_A - V_B = R I_shunt.
 */
typedef struct
{
    InternalCircuitType type;

    /*
     * Узлы открытой Stage-2 сетки.
     */
    size_t node_a;
    size_t node_b;

    /*
     * Сопротивление шунта [Ohm].
     *
     * R = 0:
     *
     *     идеальное короткое замыкание.
     */
    double resistance;


    /*
     * --------------------------------------------------------
     * ЗАПРОШЕННЫЕ ФИЗИЧЕСКИЕ КООРДИНАТЫ
     * --------------------------------------------------------
     *
     * Это координаты, заданные пользователем.
     *
     * Они особенно полезны в двухэтапной задаче:
     *
     * закрытая и открытая сетки могут иметь разные
     * номера узлов.
     */
    double requested_a[3];
    double requested_b[3];


    /*
     * --------------------------------------------------------
     * ФАКТИЧЕСКИЕ КООРДИНАТЫ УЗЛОВ СЕТКИ
     * --------------------------------------------------------
     */
    double mapped_a[3];
    double mapped_b[3];


    /*
     * Расстояние от заданной точки до выбранного узла [m].
     */
    double mapping_error_a;
    double mapping_error_b;


    /*
     * Флаги.
     */
    int nodes_mapped;
    int initialized;

} InternalCircuitShunt;


/*
 * ============================================================
 * DEFAULT
 * ============================================================
 *
 * По умолчанию:
 *
 *     R = 1000 Ohm.
 *
 * Это то же стартовое значение, которое использовалось
 * в Python Stage 2.
 */
void internal_circuit_shunt_default(
    InternalCircuitShunt *shunt
);


/*
 * ============================================================
 * ЗАДАНИЕ УЗЛОВ НАПРЯМУЮ
 * ============================================================
 *
 * Используется, если номера узлов уже известны.
 */
int internal_circuit_shunt_set_nodes(
    InternalCircuitShunt *shunt,
    const Mesh *mesh,
    size_t node_a,
    size_t node_b
);


/*
 * ============================================================
 * ЗАДАНИЕ ПО ФИЗИЧЕСКИМ КООРДИНАТАМ
 * ============================================================
 *
 * Для каждой заданной точки ищется ближайший узел сетки.
 *
 * Это предпочтительный вариант для нашей двухэтапной задачи:
 *
 *     Stage 1 mesh != Stage 2 mesh.
 *
 *
 * point_a, point_b:
 *
 *     x y z [m].
 */
int internal_circuit_shunt_map_coordinates(
    InternalCircuitShunt *shunt,
    const Mesh *mesh,
    const double point_a[3],
    const double point_b[3]
);


/*
 * ============================================================
 * ПОИСК БЛИЖАЙШЕГО УЗЛА
 * ============================================================
 */
int internal_circuit_find_nearest_node(
    const Mesh *mesh,
    const double point[3],
    size_t *node_id,
    double *distance
);


/*
 * ============================================================
 * ПРОВЕРКА
 * ============================================================
 */
int internal_circuit_shunt_validate(
    const InternalCircuitShunt *shunt,
    size_t n_nodes
);


/*
 * ============================================================
 * ВЕКТОР b_R
 * ============================================================
 *
 * Строит:
 *
 *     b_R[A] = -1
 *     b_R[B] = +1
 *
 * остальные элементы:
 *
 *     0.
 *
 *
 * Размер:
 *
 *     Nv.
 *
 *
 * Память для b_R должна быть выделена вызывающей стороной:
 *
 *     double *b_R = malloc(Nv * sizeof(double));
 */
int internal_circuit_build_incidence_vector(
    const InternalCircuitShunt *shunt,
    size_t n_nodes,
    double *b_R
);


/*
 * ============================================================
 * НАПРЯЖЕНИЕ НА ШУНТЕ
 * ============================================================
 *
 * Определение:
 *
 *     V_shunt = V_A - V_B.
 *
 *
 * Это именно то определение, которое использовалось
 * в Python-коде.
 */
double internal_circuit_shunt_voltage(
    const InternalCircuitShunt *shunt,
    const double *node_voltage,
    size_t n_nodes
);


/*
 * ============================================================
 * ЗАКОН ОМА
 * ============================================================
 *
 * Для R > 0:
 *
 *     I_shunt =
 *
 *     (V_A - V_B) / R.
 *
 *
 * Эта функция нужна только для диагностики.
 *
 * В основной Stage-2 PEEC-системе I_shunt будет
 * ОТДЕЛЬНОЙ НЕИЗВЕСТНОЙ, а не вычисляться этой функцией.
 *
 *
 * Для R = 0 функция возвращает NAN, потому что у идеального
 * КЗ ток определяется системой, а не законом V/R.
 */
double internal_circuit_shunt_current_from_voltage(
    const InternalCircuitShunt *shunt,
    const double *node_voltage,
    size_t n_nodes
);


/*
 * ============================================================
 * НЕВЯЗКА УРАВНЕНИЯ ШУНТА
 * ============================================================
 *
 * Проверяем:
 *
 *     V_A - V_B - R I_shunt = 0.
 *
 *
 * Возвращаемая величина:
 *
 *     [V].
 */
double internal_circuit_shunt_residual(
    const InternalCircuitShunt *shunt,
    const double *node_voltage,
    size_t n_nodes,
    double shunt_current
);


/*
 * ============================================================
 * ПЕЧАТЬ
 * ============================================================
 */
void internal_circuit_shunt_print_info(
    const InternalCircuitShunt *shunt
);


#endif
