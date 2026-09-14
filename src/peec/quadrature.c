#include "quadrature.h"

#include <math.h>
#include <stdio.h>


/*
 * ============================================================
 * ВСПОМОГАТЕЛЬНАЯ ГЕОМЕТРИЯ
 * ============================================================
 */


/*
 * Евклидово расстояние между двумя точками в R^3.
 */
static double distance3(const double a[3], const double b[3])
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];

    return sqrt(dx * dx + dy * dy + dz * dz);
}


/*
 * Векторное произведение:
 *
 *     result = a x b.
 */
static void cross3(const double a[3], const double b[3], double result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}


/*
 * Евклидова норма трехмерного вектора.
 */
static double norm3(const double a[3])
{
    return sqrt(
        a[0] * a[0]
        + a[1] * a[1]
        + a[2] * a[2]
    );
}


/*
 * Извлекает координаты вершины k из QuadPatch.
 *
 * В QuadPatch координаты хранятся как:
 *
 *     x0 y0 z0 x1 y1 z1 ...
 */
static void patch_get_point(
    const QuadPatch *patch,
    int k,
    double p[3]
)
{
    p[0] = patch->points[3 * k + 0];
    p[1] = patch->points[3 * k + 1];
    p[2] = patch->points[3 * k + 2];
}


/*
 * Геометрическое сравнение двух точек.
 *
 * Небольшой tolerance нужен для устойчивости при чтении
 * сеток из различных форматов.
 */
static int points_equal(const double a[3], const double b[3])
{
    const double tolerance = 1.0e-12;

    return distance3(a, b) <= tolerance;
}


/*
 * ============================================================
 * GAUSS-LEGENDRE
 * ============================================================
 *
 * Коэффициенты записаны явно, чтобы проект не зависел
 * от внешних библиотек.
 */
int gauss_legendre_rule(int order, double *x, double *w)
{
    if (x == NULL || w == NULL)
    {
        return -1;
    }


    if (order == 2)
    {
        x[0] = -0.57735026918962576451;
        x[1] =  0.57735026918962576451;

        w[0] = 1.0;
        w[1] = 1.0;

        return 0;
    }


    if (order == 3)
    {
        x[0] = -0.77459666924148337704;
        x[1] =  0.0;
        x[2] =  0.77459666924148337704;

        w[0] = 0.55555555555555555556;
        w[1] = 0.88888888888888888889;
        w[2] = 0.55555555555555555556;

        return 0;
    }


    if (order == 4)
    {
        x[0] = -0.86113631159405257522;
        x[1] = -0.33998104358485626480;
        x[2] =  0.33998104358485626480;
        x[3] =  0.86113631159405257522;

        w[0] = 0.34785484513745385737;
        w[1] = 0.65214515486254614263;
        w[2] = 0.65214515486254614263;
        w[3] = 0.34785484513745385737;

        return 0;
    }


    if (order == 5)
    {
        x[0] = -0.90617984593866399280;
        x[1] = -0.53846931010568309104;
        x[2] =  0.0;
        x[3] =  0.53846931010568309104;
        x[4] =  0.90617984593866399280;

        w[0] = 0.23692688505618908751;
        w[1] = 0.47862867049936646804;
        w[2] = 0.56888888888888888889;
        w[3] = 0.47862867049936646804;
        w[4] = 0.23692688505618908751;

        return 0;
    }


    if (order == 6)
    {
        x[0] = -0.93246951420315202781;
        x[1] = -0.66120938646626451366;
        x[2] = -0.23861918608319690863;
        x[3] =  0.23861918608319690863;
        x[4] =  0.66120938646626451366;
        x[5] =  0.93246951420315202781;

        w[0] = 0.17132449237917034504;
        w[1] = 0.36076157304813860757;
        w[2] = 0.46791393457269104739;
        w[3] = 0.46791393457269104739;
        w[4] = 0.36076157304813860757;
        w[5] = 0.17132449237917034504;

        return 0;
    }


    if (order == 8)
    {
        x[0] = -0.96028985649753623168;
        x[1] = -0.79666647741362673959;
        x[2] = -0.52553240991632898582;
        x[3] = -0.18343464249564980494;
        x[4] =  0.18343464249564980494;
        x[5] =  0.52553240991632898582;
        x[6] =  0.79666647741362673959;
        x[7] =  0.96028985649753623168;

        w[0] = 0.10122853629037625915;
        w[1] = 0.22238103445337447054;
        w[2] = 0.31370664587788728734;
        w[3] = 0.36268378337836198297;
        w[4] = 0.36268378337836198297;
        w[5] = 0.31370664587788728734;
        w[6] = 0.22238103445337447054;
        w[7] = 0.10122853629037625915;

        return 0;
    }


    if (order == 10)
    {
        x[0] = -0.97390652851717172008;
        x[1] = -0.86506336668898451073;
        x[2] = -0.67940956829902440623;
        x[3] = -0.43339539412924719080;
        x[4] = -0.14887433898163121088;
        x[5] =  0.14887433898163121088;
        x[6] =  0.43339539412924719080;
        x[7] =  0.67940956829902440623;
        x[8] =  0.86506336668898451073;
        x[9] =  0.97390652851717172008;

        w[0] = 0.06667134430868813759;
        w[1] = 0.14945134915058059315;
        w[2] = 0.21908636251598204399;
        w[3] = 0.26926671930999635509;
        w[4] = 0.29552422471475287017;
        w[5] = 0.29552422471475287017;
        w[6] = 0.26926671930999635509;
        w[7] = 0.21908636251598204399;
        w[8] = 0.14945134915058059315;
        w[9] = 0.06667134430868813759;

        return 0;
    }


    if (order == 12)
    {
        x[0]  = -0.98156063424671925069;
        x[1]  = -0.90411725637047485668;
        x[2]  = -0.76990267419430468704;
        x[3]  = -0.58731795428661744730;
        x[4]  = -0.36783149899818019375;
        x[5]  = -0.12523340851146891547;
        x[6]  =  0.12523340851146891547;
        x[7]  =  0.36783149899818019375;
        x[8]  =  0.58731795428661744730;
        x[9]  =  0.76990267419430468704;
        x[10] =  0.90411725637047485668;
        x[11] =  0.98156063424671925069;

        w[0]  = 0.04717533638651182719;
        w[1]  = 0.10693932599531843096;
        w[2]  = 0.16007832854334622633;
        w[3]  = 0.20316742672306592175;
        w[4]  = 0.23349253653835480876;
        w[5]  = 0.24914704581340278500;
        w[6]  = 0.24914704581340278500;
        w[7]  = 0.23349253653835480876;
        w[8]  = 0.20316742672306592175;
        w[9]  = 0.16007832854334622633;
        w[10] = 0.10693932599531843096;
        w[11] = 0.04717533638651182719;

        return 0;
    }


    if (order == 16)
    {
        x[0]  = -0.98940093499164993260;
        x[1]  = -0.94457502307323257608;
        x[2]  = -0.86563120238783174388;
        x[3]  = -0.75540440835500303390;
        x[4]  = -0.61787624440264374845;
        x[5]  = -0.45801677765722738634;
        x[6]  = -0.28160355077925891323;
        x[7]  = -0.09501250983763744019;
        x[8]  =  0.09501250983763744019;
        x[9]  =  0.28160355077925891323;
        x[10] =  0.45801677765722738634;
        x[11] =  0.61787624440264374845;
        x[12] =  0.75540440835500303390;
        x[13] =  0.86563120238783174388;
        x[14] =  0.94457502307323257608;
        x[15] =  0.98940093499164993260;

        w[0]  = 0.02715245941175409485;
        w[1]  = 0.06225352393864789286;
        w[2]  = 0.09515851168249278481;
        w[3]  = 0.12462897125553387205;
        w[4]  = 0.14959598881657673208;
        w[5]  = 0.16915651939500253819;
        w[6]  = 0.18260341504492358887;
        w[7]  = 0.18945061045506849629;
        w[8]  = 0.18945061045506849629;
        w[9]  = 0.18260341504492358887;
        w[10] = 0.16915651939500253819;
        w[11] = 0.14959598881657673208;
        w[12] = 0.12462897125553387205;
        w[13] = 0.09515851168249278481;
        w[14] = 0.06225352393864789286;
        w[15] = 0.02715245941175409485;

        return 0;
    }


    fprintf(
        stderr,
        "ERROR: unsupported Gauss-Legendre order: %d.\n",
        order
    );

    return -1;
}


/*
 * ============================================================
 * ПЛОЩАДЬ ТРЕУГОЛЬНИКА
 * ============================================================
 */
double triangle_area(
    const double p0[3],
    const double p1[3],
    const double p2[3]
)
{
    double a[3] = {
        p1[0] - p0[0],
        p1[1] - p0[1],
        p1[2] - p0[2]
    };

    double b[3] = {
        p2[0] - p0[0],
        p2[1] - p0[1],
        p2[2] - p0[2]
    };

    double c[3];

    cross3(a, b, c);

    return 0.5 * norm3(c);
}


/*
 * ============================================================
 * ПЛОЩАДЬ QUADPATCH
 * ============================================================
 */
double quad_patch_area(const QuadPatch *patch)
{
    if (patch == NULL)
    {
        return 0.0;
    }

    double p[4][3];

    for (int i = 0; i < 4; ++i)
    {
        patch_get_point(patch, i, p[i]);
    }

    /*
     * Разбиение:
     *
     *     (0,1,2)
     *     (0,2,3).
     */
    double area0 = triangle_area(p[0], p[1], p[2]);
    double area1 = triangle_area(p[0], p[2], p[3]);

    return area0 + area1;
}


/*
 * ============================================================
 * ХАРАКТЕРНЫЙ РАЗМЕР QUADPATCH
 * ============================================================
 */
double quad_patch_size(const QuadPatch *patch)
{
    if (patch == NULL)
    {
        return 0.0;
    }

    double p[4][3];

    for (int i = 0; i < 4; ++i)
    {
        patch_get_point(patch, i, p[i]);
    }

    double max_length = 0.0;

    /*
     * Проверяем только стороны quadrangle.
     */
    for (int i = 0; i < 4; ++i)
    {
        int j = (i + 1) % 4;

        double length = distance3(p[i], p[j]);

        if (length > max_length)
        {
            max_length = length;
        }
    }

    return max_length;
}


/*
 * ============================================================
 * МИНИМАЛЬНОЕ РАССТОЯНИЕ МЕЖДУ ВЕРШИНАМИ
 * ============================================================
 */
double quad_patch_min_vertex_distance(
    const QuadPatch *a,
    const QuadPatch *b
)
{
    if (a == NULL || b == NULL)
    {
        return 0.0;
    }

    double min_distance = HUGE_VAL;

    for (int i = 0; i < 4; ++i)
    {
        double pa[3];

        patch_get_point(a, i, pa);

        for (int j = 0; j < 4; ++j)
        {
            double pb[3];

            patch_get_point(b, j, pb);

            double distance = distance3(pa, pb);

            if (distance < min_distance)
            {
                min_distance = distance;
            }
        }
    }

    return min_distance;
}


/*
 * ============================================================
 * КОЛИЧЕСТВО ОБЩИХ ВЕРШИН QUADPATCH
 * ============================================================
 */
static int quad_common_vertices(
    const QuadPatch *a,
    const QuadPatch *b
)
{
    int common = 0;

    for (int i = 0; i < 4; ++i)
    {
        double pa[3];

        patch_get_point(a, i, pa);

        for (int j = 0; j < 4; ++j)
        {
            double pb[3];

            patch_get_point(b, j, pb);

            if (points_equal(pa, pb))
            {
                ++common;
                break;
            }
        }
    }

    return common;
}


/*
 * Проверяет полное геометрическое совпадение двух QuadPatch
 * без привязки к тому, являются ли они одним объектом.
 */
static int quad_patches_equal(
    const QuadPatch *a,
    const QuadPatch *b
)
{
    int matched[4] = {0, 0, 0, 0};

    for (int i = 0; i < 4; ++i)
    {
        double pa[3];

        patch_get_point(a, i, pa);

        int found = 0;

        for (int j = 0; j < 4; ++j)
        {
            if (matched[j])
            {
                continue;
            }

            double pb[3];

            patch_get_point(b, j, pb);

            if (points_equal(pa, pb))
            {
                matched[j] = 1;
                found = 1;
                break;
            }
        }

        if (!found)
        {
            return 0;
        }
    }

    return 1;
}


/*
 * ============================================================
 * КЛАССИФИКАЦИЯ QUADPATCH
 * ============================================================
 */
QuadPairType quad_patch_classify_pair(
    const QuadPatch *a,
    const QuadPatch *b,
    double near_factor
)
{
    if (a == NULL || b == NULL)
    {
        return QUAD_PAIR_FAR;
    }

    /*
     * Один физический объект.
     */
    if (a == b)
    {
        return QUAD_PAIR_SELF;
    }

    /*
     * Геометрически совпадающие области.
     */
    if (quad_patches_equal(a, b))
    {
        return QUAD_PAIR_SELF;
    }

    /*
     * Общая вершина или ребро.
     */
    if (quad_common_vertices(a, b) > 0)
    {
        return QUAD_PAIR_TOUCHING;
    }

    double distance2 = 0.0;
    for (int k = 0; k < 3; ++k) {
        double amin=a->points[k], amax=amin, bmin=b->points[k], bmax=bmin;
        for (int i=1; i<4; ++i) {
            amin=fmin(amin,a->points[3*i+k]); amax=fmax(amax,a->points[3*i+k]);
            bmin=fmin(bmin,b->points[3*i+k]); bmax=fmax(bmax,b->points[3*i+k]);
        }
        double gap=fmax(0.0,fmax(amin-bmax,bmin-amax));
        distance2 += gap*gap;
    }
    double d_min = sqrt(distance2);
    double h_a = quad_patch_size(a);
    double h_b = quad_patch_size(b);

    double h = (h_a > h_b) ? h_a : h_b;

    if (h > 0.0 && d_min < near_factor * h)
    {
        return QUAD_PAIR_NEAR;
    }

    return QUAD_PAIR_FAR;
}


/*
 * ============================================================
 * ПАРАМЕТРИЗАЦИЯ ТРЕУГОЛЬНИКА
 * ============================================================
 *
 * Используем преобразование квадрата:
 *
 *     u,v in [0,1]
 *
 * в треугольник:
 *
 *     xi  = u
 *     eta = (1-u)v.
 *
 * Точка:
 *
 *     r = p0 + xi(p1-p0) + eta(p2-p0).
 *
 * Поверхностный якобиан:
 *
 *     J = |(p1-p0)x(p2-p0)| (1-u).
 */
static void triangle_map(
    const double p0[3],
    const double p1[3],
    const double p2[3],
    double u,
    double v,
    double r[3],
    double *jacobian
)
{
    double a[3] = {
        p1[0] - p0[0],
        p1[1] - p0[1],
        p1[2] - p0[2]
    };

    double b[3] = {
        p2[0] - p0[0],
        p2[1] - p0[1],
        p2[2] - p0[2]
    };

    double xi = u;
    double eta = (1.0 - u) * v;

    r[0] = p0[0] + xi * a[0] + eta * b[0];
    r[1] = p0[1] + xi * a[1] + eta * b[1];
    r[2] = p0[2] + xi * a[2] + eta * b[2];

    double c[3];

    cross3(a, b, c);

    *jacobian = norm3(c) * (1.0 - u);
}


/*
 * ============================================================
 * ОБЫЧНЫЙ TRIANGLE-TRIANGLE ИНТЕГРАЛ
 * ============================================================
 *
 * Используется для несамосовпадающих пар.
 *
 * Для TOUCHING значение R становится малым около общей
 * границы, поэтому там используется более высокий порядок.
 */
static double triangle_pair_regular(
    const double a0[3],
    const double a1[3],
    const double a2[3],
    const double b0[3],
    const double b1[3],
    const double b2[3],
    int order
)
{
    if (order <= 0 || order > 16)
    {
        return NAN;
    }

    double x[16];
    double w[16];

    if (gauss_legendre_rule(order, x, w) != 0)
    {
        return NAN;
    }

    double integral = 0.0;

    /*
     * Четыре параметрические переменные:
     *
     *     u_a, v_a
     *     u_b, v_b.
     */
    for (int iu = 0; iu < order; ++iu)
    {
        double ua = 0.5 * (x[iu] + 1.0);
        double wua = 0.5 * w[iu];

        for (int iv = 0; iv < order; ++iv)
        {
            double va = 0.5 * (x[iv] + 1.0);
            double wva = 0.5 * w[iv];

            double ra[3];
            double jac_a = 0.0;

            triangle_map(
                a0,
                a1,
                a2,
                ua,
                va,
                ra,
                &jac_a
            );

            for (int ju = 0; ju < order; ++ju)
            {
                double ub = 0.5 * (x[ju] + 1.0);
                double wub = 0.5 * w[ju];

                for (int jv = 0; jv < order; ++jv)
                {
                    double vb = 0.5 * (x[jv] + 1.0);
                    double wvb = 0.5 * w[jv];

                    double rb[3];
                    double jac_b = 0.0;

                    triangle_map(
                        b0,
                        b1,
                        b2,
                        ub,
                        vb,
                        rb,
                        &jac_b
                    );

                    double R = distance3(ra, rb);

                    /*
                     * Для обычной Gauss-Legendre квадратуры
                     * точки не лежат на концах интервала,
                     * поэтому для корректной геометрии R > 0.
                     */
                    if (R <= 0.0)
                    {
                        continue;
                    }

                    double weight =
                    wua
                    * wva
                    * wub
                    * wvb
                    * jac_a
                    * jac_b;

                    integral += weight / R;
                }
            }
        }
    }

    return integral;
}


/*
 * ============================================================
 * DUFFY: СИНГУЛЯРНАЯ ВЕРШИНА ТРЕУГОЛЬНИКА
 * ============================================================
 *
 * Вычисляет:
 *
 *                  /
 *                 |       1
 *     F(r0) =     |    ------- dS
 *                 |    |r-r0|
 *                /T
 *
 * для:
 *
 *     T = (r0,p1,p2).
 *
 *
 * Параметризация:
 *
 *     r = r0 + u[(1-v)e1 + v e2],
 *
 * где:
 *
 *     e1 = p1-r0
 *     e2 = p2-r0.
 *
 *
 * Тогда:
 *
 *     R = u |(1-v)e1 + v e2|
 *
 * и:
 *
 *     dS = u |e1 x e2| du dv.
 *
 *
 * Таким образом:
 *
 *     dS/R =
 *
 *     |e1 x e2|
 *     ------------
 *     |(1-v)e1+v e2|
 *
 *     du dv.
 *
 *
 * Сингулярный множитель u полностью сокращается.
 */
/* Exact uniform-triangle potential, followed by adaptive outer integration.
 * Divergence theorem: integral_T 1/R = sum_edges d_e integral_e 1/R - |h Omega|.
 * Edge integrals use asinh; the solid angle uses the atan2 formula.
 * See van Oosterom (2011), doi:10.1007/s11517-011-0837-9, uniform monolayer.
 * No sampling of the singular 1/R kernel is needed for touching/self pairs.
 */
typedef struct {
    double p[3][3], unit[3][3], outward[3][3], length[3], normal[3];
    double twice_area;
} TrianglePotential;

static double dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void triangle_potential_prepare(double p[3][3], TrianglePotential *t)
{
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) {
            t->p[i][k] = p[i][k];
            t->unit[i][k] = p[(i+1)%3][k] - p[i][k];
        }
    cross3(t->unit[0], t->unit[1], t->normal);
    t->twice_area = norm3(t->normal);
    if (t->twice_area == 0.0) return;
    for (int k = 0; k < 3; ++k) t->normal[k] /= t->twice_area;
    for (int i = 0; i < 3; ++i) {
        t->length[i] = norm3(t->unit[i]);
        for (int k = 0; k < 3; ++k) t->unit[i][k] /= t->length[i];
        cross3(t->unit[i], t->normal, t->outward[i]);
    }
}

static double triangle_potential_at(const TrianglePotential *t, const double r[3])
{
    double v[3][3], radius[3];
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) v[i][k] = t->p[i][k] - r[k];
        radius[i] = norm3(v[i]);
    }
    double height = dot3(v[0], t->normal);
    double value = 0.0;
    for (int i = 0; i < 3; ++i) {
        double d = dot3(v[i], t->outward[i]);
        if (d == 0.0) continue; /* finite limit d log(d) = 0 */
        double along = dot3(v[i], t->unit[i]);
        double perpendicular = hypot(d, height);
        double edge_integral = asinh((along + t->length[i])/perpendicular)
                             - asinh(along/perpendicular);
        value += d * edge_integral;
    }
    if (height != 0.0) {
        double cross[3];
        cross3(v[1], v[2], cross);
        double denominator = radius[0]*radius[1]*radius[2]
            + dot3(v[0],v[1])*radius[2] + dot3(v[1],v[2])*radius[0]
            + dot3(v[2],v[0])*radius[1];
        double omega = 2.0 * atan2(dot3(v[0],cross), denominator);
        value -= fabs(height * omega);
    }
    return value;
}

static double triangle_potential_rule(double a[3][3], const TrianglePotential *b,
                                      int order)
{
    double x[16], w[16];
    if (gauss_legendre_rule(order,x,w) != 0) return NAN;
    double integral = 0.0;
    for (int i = 0; i < order; ++i)
        for (int j = 0; j < order; ++j) {
            double r[3], jac;
            triangle_map(a[0],a[1],a[2],0.5*(x[i]+1),0.5*(x[j]+1),r,&jac);
            integral += 0.25*w[i]*w[j]*jac*triangle_potential_at(b,r);
        }
    return integral;
}

static double triangle_potential_adaptive(double a[3][3], const TrianglePotential *b,
                                          int order, double abs_tol, int depth)
{
    double low = triangle_potential_rule(a,b,order <= 6 ? 3 : 4);
    double high = triangle_potential_rule(a,b,order);
    if (!isfinite(low) || !isfinite(high)) return NAN;
    if (fabs(high-low) <= abs_tol + 2e-5*fabs(high)) return high;
    if (depth >= 24) return NAN; /* never silently accept unconverged quadrature */
    int edge = 0;
    double longest = 0.0;
    for (int i=0; i<3; ++i) {
        double length = distance3(a[i],a[(i+1)%3]);
        if (length > longest) { longest = length; edge = i; }
    }
    int j = (edge+1)%3, k = (edge+2)%3;
    double left[3][3], right[3][3];
    for (int d=0; d<3; ++d) {
        left[0][d] = a[edge][d]; right[1][d] = a[j][d];
        left[1][d] = right[0][d] = 0.5*(a[edge][d]+a[j][d]);
        left[2][d] = right[2][d] = a[k][d];
    }
    return triangle_potential_adaptive(left,b,order,abs_tol*0.5,depth+1)
         + triangle_potential_adaptive(right,b,order,abs_tol*0.5,depth+1);
}

static double triangle_pair_kernel(double a[3][3], double b[3][3],
    int far_order, int near_order, int touching_order, int self_order, int near_pair)
{
    if (!near_pair)
        return triangle_pair_regular(a[0],a[1],a[2],b[0],b[1],b[2],far_order);
    TrianglePotential ta, tb;
    triangle_potential_prepare(a,&ta);
    triangle_potential_prepare(b,&tb);
    if (ta.twice_area == 0.0 || tb.twice_area == 0.0) return 0.0;
    int equal = 1;
    for (int i=0; i<3; ++i) {
        int found = 0;
        for (int j=0; j<3; ++j) if (points_equal(a[i],b[j])) found = 1;
        if (!found) equal = 0;
    }
    if (equal) {
        /* The triangle-overlap convolution has six sectors. In each,
         * overlap = area*(1-t)^2 and the radial integral is 1/3.
         * Hence K(T,T) = (2*area/3) sum_vertices integral_T 1/|v-r|.
         * This evaluates the entire self term analytically at any aspect ratio. */
        double sum = 0.0;
        for (int i=0; i<3; ++i) sum += triangle_potential_at(&tb,a[i]);
        return ta.twice_area * sum / 3.0;
    }
    /* The requested order is the outer high rule; adapt geometry as needed. */
    int order = near_pair == QUAD_PAIR_SELF+1 ? self_order
              : near_pair == QUAD_PAIR_TOUCHING+1 ? touching_order : near_order;
    if (order < 6) order = 6;
    double scale = ta.twice_area*tb.twice_area /
        sqrt(fmax(ta.twice_area,tb.twice_area));
    double ab = triangle_potential_adaptive(a,&tb,order,2e-5*scale,0);
    double ba = triangle_potential_adaptive(b,&ta,order,2e-5*scale,0);
    return 0.5*(ab+ba); /* enforce reciprocity without an ordering bias */
}

double quad_patch_kernel_1_over_r(
    const QuadPatch *a,
    const QuadPatch *b,
    int far_order,
    int near_order,
    int touching_order,
    int self_order,
    double near_factor
)
{
    if (a == NULL || b == NULL)
    {
        return NAN;
    }

    if (!gauss_legendre_order_supported(far_order) ||
        !gauss_legendre_order_supported(near_order) ||
        !gauss_legendre_order_supported(touching_order) ||
        !gauss_legendre_order_supported(self_order) ||
        !isfinite(near_factor) || near_factor < 0.0) return NAN;

    double pa[4][3];
    double pb[4][3];

    for (int i = 0; i < 4; ++i)
    {
        patch_get_point(a, i, pa[i]);
        patch_get_point(b, i, pb[i]);
    }

    /*
     * Два треугольника первой области.
     */
    double triangles_a[2][3][3] = {
        {
            {pa[0][0], pa[0][1], pa[0][2]},
            {pa[1][0], pa[1][1], pa[1][2]},
            {pa[2][0], pa[2][1], pa[2][2]}
        },
        {
            {pa[0][0], pa[0][1], pa[0][2]},
            {pa[2][0], pa[2][1], pa[2][2]},
            {pa[3][0], pa[3][1], pa[3][2]}
        }
    };

    /*
     * Два треугольника второй области.
     */
    double triangles_b[2][3][3] = {
        {
            {pb[0][0], pb[0][1], pb[0][2]},
            {pb[1][0], pb[1][1], pb[1][2]},
            {pb[2][0], pb[2][1], pb[2][2]}
        },
        {
            {pb[0][0], pb[0][1], pb[0][2]},
            {pb[2][0], pb[2][1], pb[2][2]},
            {pb[3][0], pb[3][1], pb[3][2]}
        }
    };

    QuadPairType pair_type = quad_patch_classify_pair(
        a,
        b,
        near_factor
    );

    int force_near = pair_type == QUAD_PAIR_FAR ? 0 : (int)pair_type + 1;

    double integral = 0.0;

    /*
     * Полный quad-quad интеграл:
     *
     *     T0a x T0b
     *     T0a x T1b
     *     T1a x T0b
     *     T1a x T1b.
     */
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            double value = triangle_pair_kernel(
                triangles_a[i],
                triangles_b[j],
                far_order,
                near_order,
                touching_order,
                self_order,
                force_near
            );

            if (!isfinite(value))
            {
                return NAN;
            }

            integral += value;
        }
    }

    return integral;
}
