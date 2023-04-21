#include <valarray>
#include <functional>
#include <stdexcept>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <cuba.h>
#include "point.cpp"
using namespace std;

// CONSTANTS
double mu_0 = 4*M_PI*1e-7;

//HELPFUL CLASSES
/*
 * This section of code includes methods for evaluating 2D and 3D integrals. The code is partially adapted from
 * https://stackoverflow.com/questions/40157295/efficient-double-integration-in-c-with-gsl.
 */
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// Simple RAII wrapper
class IntegrationWorkspace {
    gsl_integration_workspace * wsp;
public:
    IntegrationWorkspace(const size_t n=3*pow(10,8)):
            wsp(gsl_integration_workspace_alloc(n)) {}
    ~IntegrationWorkspace() { gsl_integration_workspace_free(wsp); }
    operator gsl_integration_workspace*() { return wsp;}
};
// Build gsl_function from lambda
template <typename F>
class gsl_function_pp: public gsl_function {
    const F func;
    static double invoke(double x, void *params) {
        return static_cast<gsl_function_pp*>(params)->func(x);
    }
public:
    gsl_function_pp(const F& f) : func(f) {
        function = &gsl_function_pp::invoke; //inherited from gsl_function
        params   = this;                     //inherited from gsl_function
    }
    operator gsl_function*(){return this;}
};
// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func) {
    return gsl_function_pp<F>(func);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// HI-FI ///////////////////////////////////////////////////////////

struct BData{
    Point r = Point(); //the Point at which we wish to determine the magnetic field
    int axis; //axis is an integer {0,1,2} that corresponds to {x,y,z}
};

/*
 * Intgrand corresponding to the function 'b' below.
 */
static int b_integrand(const int *ndim, const cubareal xx[],
                       const int *ncomp, cubareal ff[], void *data) {
    // Note that the Cuhre integration routine integrates over a unit cube, so each component of r_prime (below)
    // must be multiplied by a certain factor to ensure the integral has the proper domain.
    BData *u = (struct BData*)data;
    int axis = u->axis;
    Wire wire = u->r.get_wire();
    Point r = u->r;
    Point r_prime(r.get_wire().get_a()*xx[0], 2*M_PI*xx[1], 2*M_PI*xx[2], wire);

    valarray<double> e1 = wire.e1(r_prime.get_phi());
    double curvature = wire.kappa(r_prime.get_phi());
    double jacobian = r_prime.get_s()*(1-curvature*r_prime.get_s()*cos(r_prime.get_theta()))
                      * wire.rc_norm_firstder(r_prime.get_phi());

    double denominator = pow(r_prime.distance(r)+1e-100,3);
    double numerator;
    if(axis == 0)       {numerator = e1[1]*r_prime.dz(r)-e1[2]*r_prime.dy(r);}
    else if(axis == 1)  {numerator = e1[2]*r_prime.dx(r)-e1[0]*r_prime.dz(r);}
    else                {numerator = e1[0]*r_prime.dy(r)-e1[1]*r_prime.dx(r);}

    ff[0] = jacobian*numerator/denominator; //integrand to return
}

/*
 * This high-fidelity method returns the value of the magnetic field at a Point r along an axis {x,y,z} corresponding to
 * an int {0,1,2}, respectively. Good choices of tolerance are 1e-5.
 */
double b(Point r, int axis, double epsrel, double epsabs){
    int nregions, neval, fail;
    cubareal integral[1], error[1], prob[1];
    struct BData data;
    data.r = r;
    data.axis = axis;

    double prefactors = r.get_wire().get_I()*mu_0/r.get_wire().get_a();
    Cuhre(3, 1, b_integrand, &data, 1, epsrel, epsabs, 0,0,
          1e8, 1,nullptr, nullptr,&nregions, &neval, &fail, integral,
          error, prob);
    double result = prefactors*integral[0];
    return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// LO-FI ////////////////////////////////////////////////////////////
/*
 * This method returns the 1D regularized Biot-Savart Law along an axis {x,y,z} corresponding to an int {0,1,2},
 * respectively. This method evaluates the 1D integral in one of three ways:
 *      (1) applying Gauss-Legendre integration to the unmodified integral,
 *      (2) applying Gauss-Legendre integration to the modified integral, or
 *      (3) using an adaptive routine to accurately evaluate the integral.
 * The method of choice can be selected with the integer 'key', where (1=unmodified, 2=modified, 3=adaptive). If the user
 * selects a Gauss-Legendre routine, they may specify the number of grid points 'n_points'. Otherwise, this variable
 * is unnecessary.
 */
double b_reg(Wire wire, double phi, int axis, int key, int n_points) {
    double result;

    if (key ==1){
        auto b_approx = make_gsl_function([&](double phi_prime){
            valarray<double> vector_difference = wire.rc(phi) - wire.rc(phi_prime);
            double numerator = cross_product(wire.rc_firstder(phi_prime), vector_difference)[axis];
            double denominator = pow(pow(vector_norm(vector_difference),2)
                    + pow(wire.get_a(),2)/sqrt(M_E),1.5);
            return numerator/denominator;
        });
        result = gsl_integration_glfixed(b_approx, phi, phi+2*M_PI,
                                         gsl_integration_glfixed_table_alloc(n_points));
    }
    else if (key == 2){
        // QUADRATURE CALCULATION:
        auto b_average = make_gsl_function([&](double phi_prime){
            valarray<double> vector_difference = wire.rc(phi) - wire.rc(phi_prime);
            double numerator = cross_product(wire.rc_firstder(phi_prime), vector_difference)[axis];
            double denominator = pow(pow(vector_norm(vector_difference),2)
                    +pow(wire.get_a(),2)*exp(-0.5),1.5);
            double iota = numerator/denominator;

            double prefactor_approx = cross_product(wire.r_c_second_der(phi), wire.rc_firstder(phi))[axis];
            double numerator_approx=1-cos(phi-phi_prime);
            double denominator_approx = pow(wire.get_a()*wire.get_a()*exp(-0.5)
                    +2*(1-cos(phi-phi_prime))*pow(vector_norm(wire.rc_firstder(phi)), 2), 1.5);
            double iota_approx=prefactor_approx*numerator_approx/denominator_approx;
            return iota+iota_approx;
        });
        double quadrature_result = gsl_integration_glfixed(b_average, phi, phi+2*M_PI,
                                                           gsl_integration_glfixed_table_alloc(n_points));

        // ANALYTIC CALCULATION:
        double analytic_result= pow(vector_norm(wire.rc_firstder(phi)), -3) *
                                cross_product(wire.r_c_second_der(phi), wire.rc_firstder(phi))[axis]
                                * (1+log(wire.get_a()/(8*exp(0.25)* vector_norm(wire.rc_firstder(phi)))));

        // FINAL RESULT:
        result = quadrature_result+analytic_result;
    }
    else if (key == 3){
        auto b_approx = make_gsl_function([&](double phi_prime){
            valarray<double> vector_difference = wire.rc(phi) - wire.rc(phi_prime);
            double numerator = cross_product(wire.rc_firstder(phi_prime), vector_difference)[axis];
            double denominator = pow(pow(vector_norm(vector_difference),2)
                    + pow(wire.get_a(),2)/sqrt(M_E),1.5);
            return numerator/denominator;
        });

        double abserr;
        gsl_integration_qags(b_approx, phi, phi+2*M_PI, 1e-8, 1e-8,
                             100,IntegrationWorkspace(100),&result, &abserr);
    }
    else{throw std::invalid_argument( "key must be in {1,2,3}");}

    double prefactors = mu_0*wire.get_I()/(4*M_PI);
    return prefactors * result;
}

/*
 * This method returns the 1D integral solution for the magnetic field at a Point 'r' if s<=a along an axis {x,y,z}
 * corresponding to an int {0,1,2}, respectively. This method automatically uses an adaptive integration routine.
 * NOTE: if s>a, the method throws an error.
 */
double b_1D(Point r, int axis){
    Wire coil = r.get_wire();
    double a = coil.get_a();
    double s = r.get_s();
    double theta = r.get_theta();
    double phi = r.get_phi();

    if(s>a){throw std::invalid_argument( "s must be <= a");}

    //Infinite cylinder solution for B
    valarray<double> b_cylinder=mu_0*coil.get_I()*s/(2*M_PI*a*a)*(-sin(theta)* coil.e2(phi) + cos(theta) * coil.e3(phi));

    //Extra terms in solution for B
    valarray<double> b_extra=-s*s/(a*a*2)*sin(2*theta)* coil.e2(phi)
            + (0.75 + s * s / (a * a) * (cos(2 * theta) / 2 - 1) * coil.e3(phi));
    b_extra*= mu_0 * coil.get_I() * coil.kappa(phi) / (8 * M_PI);

    if(axis==0){return b_cylinder[0] + b_extra[0] + b_reg(coil, phi, 0, 3, -1);}
    else if(axis == 1){return b_cylinder[1] + b_extra[1] + b_reg(coil, phi, 1, 3, -1);}
    else{return b_cylinder[2] + b_extra[2] + b_reg(coil, phi, 2, 3, -1);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// MAX|B| CALCULATIONS ////////////////////////////////////////////////////
struct ModBData{Wire w; double phi;};

/*
 * Returns -|B| in a hi-fi form compatible with GSL solvers. The quantity 'params' should refer to a struct of the type
 * ModBData.
 */
double neg_modB (const gsl_vector *v, void *params)
{
    ModBData *u = (struct ModBData*)params;
    Wire wire = u->w;
    double phi = u->phi;
    double s = gsl_vector_get(v, 0);
    double theta = gsl_vector_get(v, 1);

    Point r(s, theta, phi, wire);
    double bx = b(r, 0, 1e-5, 1e-5);
    double by = b(r, 1, 1e-5, 1e-5);
    double bz = b(r, 2, 1e-5, 1e-5);

    return -sqrt(bx*bx + by*by + bz*bz);
}

/*
 * This high-fidelity method returns the maximum value of |B| along the slice defined by the toroidal angle phi.
 * It does so by minimizing -|B| using the Nelder-Mead Simplex algorithm, where B is determined using the high-fidelity
 * method. Smaller simplexes result in tighter tolerances, a reasonable value is 0.001.
 */
double max_modB(Wire wire, double phi, double simplex_size){
    ModBData par; par.w = wire; par.phi = phi;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point, (s,theta)=(a,0) */
    x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, wire.get_a());
    gsl_vector_set (x, 1, 0);

    /* Set initial step sizes to a/2 */
    ss = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 0.5*wire.get_a());

    /* Initialize method and iterate */
    minex_func.n = 2; minex_func.f = neg_modB; minex_func.params =&par;

    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do{
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status)
        break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, simplex_size); //initially 1e-2
    }
    while (status == GSL_CONTINUE && iter < 100);

    gsl_vector_free(x); gsl_vector_free(ss); gsl_multimin_fminimizer_free (s);
    return -1*(s->fval);
}

/*
 * This low-fidelity method returns the maximum value of |B| along the slice defined by the toroidal angle phi. This
 * method is fully analytic in the sense that the magnetic field is determined using the 1D method and max|B| was
 * approximated from this analytically.
 */
double max_modB_fullyanalytic(Wire wire, double phi){
    // Magnetic field in cartesian coordinates
    double bx = b_reg(wire, phi, 0, 3, -1);
    double by = b_reg(wire, phi, 1, 3, -1);
    double bz = b_reg(wire, phi, 2, 3, -1);
    valarray<double> b_cartesian = {bx, by, bz};

    //Magnetic field in frenet-serret coordinates
    double b1 = dot_product(b_cartesian, wire.e1(phi));
    double b2 = dot_product(b_cartesian, wire.e2(phi));
    double b3 = dot_product(b_cartesian, wire.e3(phi));
    double b_par = sqrt(b2*b2 + b3*b3); //definition of b_parallel

    //finding the maximum magnetic field
    double term1 = b1;
    double term2 = b2+mu_0*wire.get_I()/(2*M_PI*wire.get_a())*(b2/b_par+
            wire.kappa(phi) * wire.get_a() * b2 * b3 / (8 * b_par * b_par));
    double term3 = b3+mu_0*wire.get_I()/(2*M_PI*wire.get_a())*(b3/b_par+
            wire.kappa(phi) * wire.get_a() / 8 * (-0.5 + (b3 * b3 - b2 * b2) / (b_par * b_par)));

    return sqrt(term1*term1 + term2*term2 + term3*term3);
}

/*
 * Returns -|B| in a lo-fi form compatible with GSL solvers. The quantity 'params' should refer to a struct of the type
 * ModBData.
 */
double neg_modB_1D (const gsl_vector *v, void *params)
{
    ModBData *u = (struct ModBData*)params;
    Wire wire = u->w;
    double s = gsl_vector_get(v, 0);
    if (s>wire.get_a()){s=wire.get_a();} // b_1D doesn't take s>a
    double theta = gsl_vector_get(v, 1);
    double phi = u->phi;

    Point r(s, theta, phi, wire);
    double bx = b_1D(r, 0);
    double by = b_1D(r, 1);
    double bz = b_1D(r, 2);

    return -sqrt(bx*bx + by*by + bz*bz);
}

/*
 * This low-fidelity method returns the maximum value of |B| along the slice defined by the toroidal angle phi. This
 * method is semi-analytic in the sense that the magnetic field is determined using the 1D method and max|B| was
 * determined with the Nelder-Mead Simplex minimization routine. The default desired size of the simplex is set to be
 * 0.001 and gives reasonable results; smaller simplexes result in tighter tolerances.
 */
double max_modB_semianalytic(Wire wire, double phi, double simplex_size){
    ModBData par; par.w = wire; par.phi = phi;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point, (s,theta)=(a,0) */
    x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, wire.get_a());
    gsl_vector_set (x, 1, 0);

    /* Set initial step sizes to a/2 */
    ss = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 0.5*wire.get_a());

    /* Initialize method and iterate */
    minex_func.n = 2; minex_func.f = neg_modB_1D; minex_func.params =&par;

    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do{
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, simplex_size); //initially 1e-2
    }
    while (status == GSL_CONTINUE && iter < 100);

    gsl_vector_free(x); gsl_vector_free(ss); gsl_multimin_fminimizer_free (s);
    return -1*(s->fval);
}