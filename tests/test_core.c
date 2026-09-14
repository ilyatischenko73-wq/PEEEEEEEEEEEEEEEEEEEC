#include "aperture_coupling.h"
#include "complex_lu.h"
#include "potential.h"
#include "harmonic.h"
#include "quadrature.h"
#include "transient.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Acceptance targets for these normalized, analytic reference problems. */
#define ALGEBRA_TOL 1.0e-11
#define QUADRATURE_REL_TOL 1.0e-3

static int test_real_lu(void)
{
    const double a[] = {0, 2, 1, 1, -2, -3, 2, 3, 1};
    const double exact[] = {1, -2, 0.5};
    double rhs[3] = {0}, x[3] = {0};
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j) rhs[i] += a[3*i+j]*exact[j];
    LUOptions options;
    LUFactorization lu;
    lu_options_default(&options, 1);
    lu_init(&lu);
    int status = lu_factorize(a, 3, &options, &lu);
    if (status == 0) status = lu_solve(&lu, rhs, x);
    double error = 0;
    for (size_t i = 0; i < 3; ++i) {
        if (!isfinite(x[i])) status = -1;
        error = fmax(error, fabs(x[i]-exact[i]));
    }
    printf("real LU: max solution error = %.9e\n", error);
    lu_free(&lu);
    return status == 0 && error < ALGEBRA_TOL;
}

static int test_complex_lu(void)
{
    const Complex a[] = {{0,0}, {2,1}, {1,-1}, {3,0.5}};
    const Complex exact[] = {{1,2}, {-0.5,0.75}};
    Complex rhs[2] = {{0,0}, {0,0}}, x[2] = {{0,0}, {0,0}};
    /* Independent scalar arithmetic for constructing the reference RHS. */
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            rhs[i].re += a[2*i+j].re*exact[j].re-a[2*i+j].im*exact[j].im;
            rhs[i].im += a[2*i+j].re*exact[j].im+a[2*i+j].im*exact[j].re;
        }
    }
    ComplexLUOptions options;
    ComplexLUFactorization lu;
    complex_lu_options_default(&options, 1);
    complex_lu_init(&lu);
    int status = complex_lu_factorize(a, 2, &options, &lu);
    if (status == 0) status = complex_lu_solve(&lu, rhs, x);
    double error = 0;
    for (size_t i = 0; i < 2; ++i) {
        if (!isfinite(x[i].re) || !isfinite(x[i].im)) status = -1;
        error = fmax(error, hypot(x[i].re-exact[i].re, x[i].im-exact[i].im));
    }
    printf("complex LU: max solution error = %.9e\n", error);
    complex_lu_free(&lu);
    return status == 0 && error < ALGEBRA_TOL;
}

/* One branch, two nodes, nonzero mutual electrostatic coefficient.
 * L=1, R=0, P=[[1,.25],[.25,1]], A=[-1,1], F=identity.
 * I(0)=1, Vc(0)=0 -> omega=sqrt(1.5), I=cos(omega*t),
 * Vc=[-sin(omega*t)/omega, +sin(omega*t)/omega].
 */
static double oscillator_energy(double current, const double v[2])
{
    return 0.5*current*current
        + 0.5*(v[0]*v[0]+v[1]*v[1]+0.5*v[0]*v[1]);
}

static int oscillator_run(TransientScheme scheme, int steps, double *error)
{
    double ld[] = {1}, rd[] = {0}, ad[] = {-1,1};
    double pd[] = {1,0.25,0.25,1};
    InductanceMatrix l = {.n=1, .data=ld};
    ResistanceMatrix r = {.n=1, .diagonal=rd};
    IncidenceMatrix a = {.rows=1, .cols=2, .data=ad};
    PotentialMatrix p = {.n=2, .data=pd};
    TransientSystem system;
    transient_system_init(&system);
    double current[] = {1}, v[] = {0,0};
    const double dt = 1.0/steps;
    double previous_energy = 0.5, max_energy_error = 0;
    int ok = 0;
    if (transient_system_create_stage1(&system, &l, &r, &a, &p, scheme, dt)
        || transient_factorize(&system, 1)) goto cleanup;
    for (int n = 0; n < steps; ++n) {
        if (transient_build_rhs_stage1(&system, current, v, NULL,NULL,NULL,NULL)
            || transient_solve_step(&system)
            || transient_extract_stage1_solution(&system, current, v)) goto cleanup;
        const double energy = oscillator_energy(current[0], v);
        if (!isfinite(energy) || !isfinite(v[0]) || !isfinite(v[1])) goto cleanup;
        if (fabs(v[0]+v[1]) > ALGEBRA_TOL) goto cleanup; /* total charge */
        if (scheme == TRANSIENT_TRAPEZOIDAL && fabs(energy-0.5) > ALGEBRA_TOL)
            goto cleanup;
        if (scheme == TRANSIENT_BACKWARD_EULER && energy > previous_energy+ALGEBRA_TOL)
            goto cleanup;
        max_energy_error = fmax(max_energy_error, fabs(energy-0.5));
        previous_energy = energy;
    }
    const double omega = sqrt(1.5), ref_v = sin(omega)/omega;
    *error = fmax(fabs(current[0]-cos(omega)),
        fmax(fabs(v[0]+ref_v), fabs(v[1]-ref_v)));
    printf("dt=%.6g error=%.9e energy_drift=%.9e\n", dt,*error,max_energy_error);
    ok = isfinite(*error) && *error > 0;
cleanup:
    transient_system_free(&system);
    return ok;
}

static int test_time_order(TransientScheme scheme)
{
    double errors[3];
    for (int k = 0; k < 3; ++k)
        if (!oscillator_run(scheme, 10*(1<<k), &errors[k])) return 0;
    const double expected = scheme == TRANSIENT_TRAPEZOIDAL ? 4.0 : 2.0;
    for (int k = 1; k < 3; ++k) {
        double ratio = errors[k-1]/errors[k];
        printf("error ratio = %.8f; expected about %.1f\n", ratio, expected);
        if (!isfinite(ratio) || fabs(ratio-expected) > 0.1*expected) return 0;
    }
    return 1;
}

static int test_be(void) { return test_time_order(TRANSIENT_BACKWARD_EULER); }
static int test_trapezoidal(void) { return test_time_order(TRANSIENT_TRAPEZOIDAL); }

static int test_shunt_energy(void)
{
    double ld[]={1}, rd[]={0}, ad[]={-1,1}, pd[]={1,0.25,0.25,1};
    double xyz[]={0,0,0,1,0,0};
    InductanceMatrix l={.n=1,.data=ld};
    ResistanceMatrix r={.n=1,.diagonal=rd};
    IncidenceMatrix a={.rows=1,.cols=2,.data=ad};
    PotentialMatrix p={.n=2,.data=pd};
    Mesh mesh={.n_nodes=2,.xyz=xyz};
    InternalCircuitShunt shunt;
    internal_circuit_shunt_default(&shunt);
    shunt.resistance=2;
    if (internal_circuit_shunt_set_nodes(&shunt,&mesh,0,1)) return 0;
    TransientSystem system;
    transient_system_init(&system);
    const double dt=0.025;
    double current[]={0}, v[]={1,-1}, ishunt=0.75;
    double max_balance=0, max_ohm=0;
    int ok=0;
    if (transient_system_create_stage2(&system,&l,&r,&a,&p,&shunt,
            TRANSIENT_TRAPEZOIDAL,dt) || transient_factorize(&system,1)) goto cleanup;
    for (int step=0; step<40; ++step) {
        const double previous_energy=oscillator_energy(current[0],v);
        const double previous_ishunt=ishunt;
        if (transient_build_rhs_stage2(&system,current,v,ishunt,NULL,NULL,NULL,NULL)
            || transient_solve_step(&system)
            || transient_extract_stage2_solution(&system,current,v,&ishunt)) goto cleanup;
        double energy=oscillator_energy(current[0],v);
        double average_ishunt=0.5*(ishunt+previous_ishunt);
        double balance=energy-previous_energy+dt*shunt.resistance*average_ishunt*average_ishunt;
        double ohm=0.75*(v[0]-v[1])-shunt.resistance*ishunt;
        if (!isfinite(balance) || !isfinite(ohm) || !isfinite(energy)) goto cleanup;
        max_balance=fmax(max_balance,fabs(balance));
        max_ohm=fmax(max_ohm,fabs(ohm));
        if (fabs(v[0]+v[1]) > ALGEBRA_TOL) goto cleanup;
    }
    printf("shunt energy-balance error=%.9e, Ohm error=%.9e\n",max_balance,max_ohm);
    ok=max_balance<ALGEBRA_TOL && max_ohm<ALGEBRA_TOL;
cleanup:
    transient_system_free(&system);
    return ok;
}

/* Exact rectangle self integral. Using the displacement convolution:
 * G = 4 int_0^a int_0^b (a-x)(b-y)/sqrt(x*x+y*y) dy dx.
 * Long double reduces cancellation in its closed form.
 */
static long double rectangle_reference(long double a, long double b)
{
    long double h=hypotl(a,b);
    long double i00=a*asinhl(b/a)+b*asinhl(a/b);
    long double ix=(b*h+a*a*asinhl(b/a)-b*b)/2;
    long double iy=(a*h+b*b*asinhl(a/b)-a*a)/2;
    long double ixy=(h*h*h-a*a*a-b*b*b)/3;
    return 4*(a*b*i00-b*ix-a*iy+ixy);
}

static int test_rectangle(double a)
{
    QuadPatch patch={.points={0,0,0,a,0,0,a,1,0,0,1,0},.owner_id=0,.cell_id=0};
    /* Test the current production defaults, including future changes to them. */
    PotentialOptions options;
    potential_options_default(&options,1);
    double actual=quad_patch_kernel_1_over_r(&patch,&patch,
        options.far_order,options.near_order,options.touching_order,
        options.self_order,options.near_factor);
    long double expected=rectangle_reference(a,1);
    double relative=(double)fabsl(((long double)actual-expected)/expected);
    printf("aspect=%.0f exact=%.12Lg actual=%.12g relative_error=%.6g limit=%.6g\n",
        a,expected,actual,relative,QUADRATURE_REL_TOL);
    return isfinite(actual) && actual>0 && isfinite(relative) && relative<=QUADRATURE_REL_TOL;
}

static int test_square(void) { return test_rectangle(1); }
static int test_rectangle10(void) { return test_rectangle(10); }
static int test_rectangle100(void) { return test_rectangle(100); }

static int aperture_case(int correct_geometry)
{
    double closed_xyz[]={0,0,0,1,0,0};
    double open_xyz[]={0,0,0};
    if (!correct_geometry) open_xyz[0]=100;
    int edges[]={0,1};
    Mesh closed={.n_nodes=2,.n_edges=1,.xyz=closed_xyz,.edges=edges};
    Mesh open={.n_nodes=1,.xyz=open_xyz};
    ApertureNodeMapEntry map={.closed_node=0,.open_node=0};
    size_t cover=0;
    ApertureCoupling coupling;
    aperture_coupling_init(&coupling);
    int status=aperture_coupling_build(&coupling,&closed,&open,&map,1,&cover,1);
    int ok=status!=0;
    if (correct_geometry) {
        double current=2, source=0;
        ok=status==0 && aperture_coupling_apply(&coupling,&current,&source)==0
            && isfinite(source) && fabs(source+2)<ALGEBRA_TOL;
        printf("valid aperture: status=%d, source=%.9g (expected -2)\n",status,source);
    } else {
        printf("invalid aperture: status=%d (must be nonzero)\n",status);
    }
    aperture_coupling_free(&coupling);
    return ok;
}

static int test_aperture_valid(void) { return aperture_case(1); }
static int test_aperture_invalid(void) { return aperture_case(0); }

/* Geometric identities exercise rotated panels, touching pairs, disjoint
 * panels and near parallel interiors, independently of the implementation. */
static double kernel(const QuadPatch *a, const QuadPatch *b)
{
    PotentialOptions o;
    potential_options_default(&o,1);
    return quad_patch_kernel_1_over_r(a,b,o.far_order,o.near_order,
                                     o.touching_order,o.self_order,o.near_factor);
}

static int test_quadrature_geometry(void)
{
    QuadPatch a={.points={0,0,0, 2,0,0, 2,1,0, 0,1,0}};
    QuadPatch b=a;
    /* Orthonormal axes in 3D and a translation; scale gives a cubic integral. */
    const double e[3]={0.6,0.8,0}, f[3]={-0.48,0.36,0.8};
    for (int i=0;i<4;++i)
        for (int k=0;k<3;++k)
            b.points[3*i+k]=3.7*(a.points[3*i]*e[k]+a.points[3*i+1]*f[k])+k+1;
    double exact=(double)rectangle_reference(2,1)*3.7*3.7*3.7;
    double err=fabs(kernel(&b,&b)/exact-1);
    printf("rotated/translated/scaled self relative error=%.6g\n",err);
    return isfinite(err) && err<QUADRATURE_REL_TOL;
}

static int test_quadrature_partition(void)
{
    QuadPatch a={.points={0,0,0, 1,0,0, 1,1,0, 0,1,0}}, b=a;
    for (int i=0;i<4;++i) b.points[3*i]+=1;
    double ab=kernel(&a,&b), ba=kernel(&b,&a);
    double total=kernel(&a,&a)+kernel(&b,&b)+ab+ba;
    double exact=(double)rectangle_reference(2,1);
    double err=fabs(total/exact-1);
    printf("touching partition relative error=%.6g, reciprocity=%.6g\n",err,fabs(ab-ba));
    return isfinite(err) && err<QUADRATURE_REL_TOL && fabs(ab-ba)<1e-12;
}

static int test_quadrature_near_far(void)
{
    QuadPatch a={.points={0,0,0, 1,0,0, 1,1,0, 0,1,0}}, b=a;
    for (int i=0;i<4;++i) b.points[3*i+2]=1e-7;
    double near=kernel(&a,&b), self=(double)rectangle_reference(1,1);
    double near_error=fabs(near/self-1);
    for (int i=0;i<4;++i) b.points[3*i+2]=100;
    double far=kernel(&a,&b);
    /* 1/sqrt(100^2+2) <= 1/R <= 1/100 over these two unit panels. */
    int ok=isfinite(near_error) && near_error<QUADRATURE_REL_TOL
        && far>=1/sqrt(10000.0+2) && far<=0.01;
    printf("near-to-self relative difference=%.6g, far=%.12g\n",near_error,far);
    return ok;
}

static int test_harmonic_loop_reference(void)
{
    /* Three-edge cycle: a known divergence-free loop plus a gradient drive.
     * The graph Laplacian has eigenvalue 3 on the zero-sum subspace. */
    double l[9]={1e-7,0,0,0,1e-7,0,0,0,1e-7}, r[3]={0};
    double a[9]={-1,1,0,0,-1,1,1,0,-1};
    double p[9]={1e11,0,0,0,1e11,0,0,0,1e11}, pat[9];
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) pat[3*i+j]=1e11*a[3*j+i];
    InductanceMatrix L={.n=3,.data=l}; ResistanceMatrix R={.n=3,.diagonal=r};
    IncidenceMatrix A={.rows=3,.cols=3,.data=a};
    PotentialMatrix P={.n=3,.data=p}; Matrix PAT={.rows=3,.cols=3,.data=pat};
    double worst=0;
    for (int k=0;k<3;++k) {
        double frequency=k==0?1:k==1?1e3:1e6;
        double w=2*acos(-1.0)*frequency, den=3e11-w*w*1e-7;
        Complex u[3]={{1,-w*1e-7*0.003},{-1,-w*1e-7*0.003},{0,-w*1e-7*0.003}};
        HarmonicExcitation ex={.n_edges=3,.n_nodes=3,.edge_voltage=u};
        HarmonicSystem system; HarmonicSolution sol;
        harmonic_system_init(&system); harmonic_solution_init(&sol);
        int status=harmonic_build_system(&L,&R,&A,&PAT,frequency,1,&system);
        if (!status) status=harmonic_factorize(&system,1);
        if (!status) status=harmonic_solve(&system,&P,&ex,&sol,1);
        if (status) { harmonic_system_free(&system); harmonic_solution_free(&sol); return 0; }
        double e[3]={1,-1,0}, v[3]={1e11/den,-2e11/den,1e11/den};
        for (int i=0;i<3;++i) {
            worst=fmax(worst,hypot(sol.edge_current[i].re-0.003,sol.edge_current[i].im+w*e[i]/den)/0.003);
            worst=fmax(worst,hypot(sol.node_voltage[i].re-v[i],sol.node_voltage[i].im));
        }
        harmonic_system_free(&system); harmonic_solution_free(&sol);
    }
    printf("analytic loop/gradient solution max relative error=%.9g\n",worst);
    return isfinite(worst) && worst<1e-6;
}

static int test_complex_lu_singular(void)
{
    Complex a[4]={{1,1},{2,2},{2,2},{4,4}};
    ComplexLUOptions o; ComplexLUFactorization lu;
    complex_lu_options_default(&o,1); complex_lu_init(&lu);
    int status=complex_lu_factorize(a,2,&o,&lu);
    complex_lu_free(&lu);
    printf("singular matrix status=%d (must be nonzero)\n",status);
    return status!=0;
}

typedef int (*TestFunction)(void);
typedef struct { const char *name; TestFunction run; } TestCase;

int main(int argc, char **argv)
{
    const TestCase cases[]={
        {"quadrature_geometry",test_quadrature_geometry},
        {"quadrature_partition",test_quadrature_partition},
        {"quadrature_near_far",test_quadrature_near_far},
        {"harmonic_loop_reference",test_harmonic_loop_reference},
        {"complex_lu_singular",test_complex_lu_singular},
        {"real_lu",test_real_lu},
        {"complex_lu",test_complex_lu},
        {"backward_euler_order",test_be},
        {"trapezoidal_order_energy",test_trapezoidal},
        {"shunt_energy_ohm",test_shunt_energy},
        {"quadrature_square",test_square},
        {"quadrature_rectangle10",test_rectangle10},
        {"quadrature_rectangle100",test_rectangle100},
        {"aperture_valid",test_aperture_valid},
        {"aperture_reject_geometry",test_aperture_invalid},
    };
    if (argc!=2) {
        fprintf(stderr,"Usage: %s --list | TEST_NAME\n",argv[0]);
        return 2;
    }
    for (size_t i=0;i<sizeof(cases)/sizeof(cases[0]);++i) {
        if (strcmp(argv[1],"--list")==0) puts(cases[i].name);
        else if (strcmp(argv[1],cases[i].name)==0) {
            int ok=cases[i].run();
            printf("%s %s\n",ok?"PASS":"FAIL",cases[i].name);
            return ok ? EXIT_SUCCESS : EXIT_FAILURE;
        }
    }
    if (strcmp(argv[1],"--list")==0) return EXIT_SUCCESS;
    fprintf(stderr,"Unknown test: %s\n",argv[1]);
    return 2;
}
