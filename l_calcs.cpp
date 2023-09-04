#include "f_calcs.cpp"
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// HI-FI ///////////////////////////////////////////////////////////

struct AData{
    Point r = Point(); //the Point at which we wish to determine the magnetic field
    int axis; //axis is an integer {0,1,2} that corresponds to {x,y,z}
};

/*
 * Intgrand corresponding to the function 'A' below.
 */
static int A_integrand(const int *ndim, const cubareal xx[],
                       const int *ncomp, cubareal ff[], void *data) {
    // Note that the Cuhre integration routine integrates over a unit cube, so each component of r_prime (below)
    // must be multiplied by a certain factor to ensure the integral has the proper domain.
    AData *u = (struct AData*)data;
    int axis = u->axis;
    Wire wire = u->r.get_wire();
    Point r = u->r;

    double theta = r.get_theta();
    double phi = r.get_phi();
    double s_p = r.get_wire().get_a()*xx[0];
    double theta_p = 2*M_PI*xx[1]+theta; //shift domain such that singularity is on boundary
    double phi_p = 2*M_PI*xx[2]+phi; //shift domain such that singularity is on boundary

    Point r_prime(s_p, theta_p, phi_p, wire);

    valarray<double> e1 = wire.e1(phi_p);
    double curvature = wire.kappa(phi_p);
    double jacobian = s_p*(1-curvature*s_p*cos(theta_p));

    double numerator = wire.rc_norm_firstder(phi_p)*e1[axis];
    double denominator = r_prime.distance(r);

    ff[0] = jacobian*numerator/denominator; //integrand to return
}


/*
 * This high-fidelity method returns the value of the vector potential at a Point r along an axis {x,y,z} corresponding to
 * an int {0,1,2}, respectively. Good choices of tolerance are 1e-5.
 */
double A(Point r, int axis, double epsrel, double epsabs){
    int nregions, neval, fail;
    cubareal integral[1], error[1], prob[1];
    struct AData data;
    data.r = r;
    data.axis = axis;

    double prefactors = r.get_wire().get_I()*mu_0;
    Cuhre(3, 1, A_integrand, &data, 1, epsrel, epsabs, 0,0,
          100000000, 1,nullptr, nullptr,&nregions, &neval, &fail, integral,
          error, prob);
    double result = prefactors*integral[0];
    return result;
}

struct LData{
    double epsrel; double epsabs;
    Wire wire;
};

/*
 * Intgrand corresponding to the function 'L' below.
 */
static int L_integrand(const int *ndim, const cubareal xx[],
                       const int *ncomp, cubareal ff[], void *data) {
    // Note that the Cuhre integration routine integrates over a unit cube, so each component of r_prime (below)
    // must be multiplied by a certain factor to ensure the integral has the proper domain.
    LData *u = (struct LData*)data;
    Wire wire = u->wire;
    double epsrel = u->epsrel;
    double epsabs = u->epsabs;

    double s_p = wire.get_a()*xx[0];
    double theta_p = 2*M_PI*xx[1];
    double phi_p = 2*M_PI*xx[2];

    Point r_prime(s_p, theta_p, phi_p, wire);

    valarray<double> e1 = wire.e1(phi_p);
    double curvature = wire.kappa(phi_p);
    double jacobian = s_p*(1-curvature*s_p*cos(theta_p))*wire.rc_norm_firstder(phi_p);

    valarray<double> vec_pot = {A(r_prime,0,epsrel, epsabs),
                                A(r_prime,1,epsrel, epsabs),
                                A(r_prime,2,epsrel, epsabs)};
    double numerator = dot_product(e1, vec_pot);

    ff[0] = jacobian*numerator;
}

/*
 * This high-fidelity method returns the value of the vector potential. Good tolerances seem to be
 * epsrel=1e-7 and epsabs=1e-5.
 */
double L(Wire wire, double epsrel, double epsabs){
    int nregions, neval, fail;
    cubareal integral[1], error[1], prob[1];
    struct LData data;
    data.wire = wire;
    data.epsrel = epsrel;
    data.epsabs = epsabs;

    double prefactors = 4*M_PI/wire.get_I()/(wire.get_a()*wire.get_a());
    Cuhre(3, 1, L_integrand, &data, 1, epsrel, epsabs, 0,0,
          100000000, 1,nullptr, nullptr,&nregions, &neval, &fail, integral,
          error, prob);
    double result = prefactors*integral[0];
    return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// LO-FI ///////////////////////////////////////////////////////////

/*
 * Intgrand corresponding to the function 'L_2D' below.
 */
static int L_2D_integrand(const int *ndim, const cubareal xx[],
                       const int *ncomp, cubareal ff[], void *data) {
    // Note that the Cuhre integration routine integrates over a unit cube, so each component of r_prime (below)
    // must be multiplied by a certain factor to ensure the integral has the proper domain.
    LData *u = (struct LData*)data;
    Wire wire = u->wire;

    double phi_u = 2*M_PI*xx[0];
    double phi_v = 2*M_PI*xx[1]+phi_u;

    Point p_u(0,0,phi_u,wire);
    Point p_v(0,0,phi_v,wire);

    double numerator = dot_product(wire.rc_firstder(phi_u),wire.rc_firstder(phi_v));
    double denominator = sqrt(pow(p_u.distance(p_v),2)+wire.get_a()*wire.get_a()/sqrt(M_E));
    double prefactors = mu_0*M_PI;

    ff[0] = prefactors*numerator/denominator;
}

/*
 * This low-fidelity method returns the value of the inductance.
 */
double L_2D(Wire wire, double epsrel, double epsabs){
    int nregions, neval, fail;
    cubareal integral[1], error[1], prob[1];
    struct LData data;
    data.wire = wire;

    Cuhre(2, 1, L_2D_integrand, &data, 1, epsrel, epsabs, 0,0,
          100000000, 1,nullptr, nullptr,&nregions, &neval, &fail, integral,
          error, prob);
    double result = integral[0];
    return result;
}

/*
 * Returns the first term in the Fourier-series based approximation for the inductance.
 */
double L_fourier(Wire w){
    double norm_a=sqrt(w.get_xc_coeff()[1]*w.get_xc_coeff()[1]+w.get_yc_coeff()[1]*w.get_yc_coeff()[1]+w.get_zc_coeff()[1]*w.get_zc_coeff()[1]);
    double norm_b=sqrt(w.get_xs_coeff()[1]*w.get_xs_coeff()[1]+w.get_ys_coeff()[1]*w.get_ys_coeff()[1]+w.get_zs_coeff()[1]*w.get_zs_coeff()[1]);
    double R=sqrt(norm_a*norm_a+norm_b*norm_b);

    double L=0;
    for(int n=1; n<w.get_order(); n++){
        double a_n=sqrt(w.get_xc_coeff()[n]*w.get_xc_coeff()[n]+w.get_yc_coeff()[n]*w.get_yc_coeff()[n]+w.get_zc_coeff()[n]*w.get_zc_coeff()[n]);
        double b_n=sqrt(w.get_xs_coeff()[n]*w.get_xs_coeff()[n]+w.get_ys_coeff()[n]*w.get_ys_coeff()[n]+w.get_zs_coeff()[n]*w.get_zs_coeff()[n]);
        L+=(a_n*a_n+b_n*b_n)*(-sqrt(2)*log(w.get_a()*w.get_a()/(2*R*R*sqrt(M_E)))-4*sqrt(2)+5*sqrt(2)*log(2)-2.77216*log(n));
    }
    L*=mu_0/(4*sqrt(2)*R);
    return L;
}

/*
 * This method returns the 2D inductance and evaluates the integral in one of two ways:
 *    (1) applying Gauss-Legendre integration to the unmodified integral,
 *    (2) applying Gauss-Legendre integration to the modified integral, or
 * The method of choice can be selected with the integer 'key', where (1=unmodified, 2=modified).
 */
double L_2D_quadrature(Wire wire, int grid_size, int key){
    if(key==1){
        //Calculating the 2D integral by 2D gauss-legendre rule
        auto outer = make_gsl_function([&](double phi){
            auto inner = make_gsl_function([&](double phi_p){
                Point r(0, 0, phi, wire);
                Point r_p(0, 0, phi_p, wire);

                double numerator_1= dot_product(wire.rc_firstder(phi),wire.rc_firstder(phi_p));
                double denominator_1=sqrt(r.distance(r_p)*r.distance(r_p)+wire.get_a()*wire.get_a()/sqrt(M_E));
                return numerator_1/denominator_1;
            });
            double result_inner = gsl_integration_glfixed(inner, phi, phi+2*M_PI,
                                                          gsl_integration_glfixed_table_alloc(grid_size));
            return result_inner;
        });
        return mu_0/(4*M_PI)*gsl_integration_glfixed(outer, 0, 2*M_PI,
                                                          gsl_integration_glfixed_table_alloc(grid_size));
    }
    else{
        //calculating hte 1D integral
        double result_1D, abserr;
        auto L_1D = make_gsl_function([&](double phi){
            return mu_0/(4*M_PI)*wire.rc_norm_firstder(phi)*(2*log(8*wire.rc_norm_firstder(phi)/wire.get_a())+0.5);
        });
        gsl_integration_qags(L_1D, 0, 2*M_PI, 1e-8, 1e-8,
                             100,IntegrationWorkspace(100),&result_1D, &abserr);

        //calculating the 2D integral
        double result_2D;
        auto outer = make_gsl_function([&](double phi){
            auto inner = make_gsl_function([&](double phi_p){
                Point r(0, 0, phi, wire);
                Point r_p(0, 0, phi_p, wire);

                //integral 1
                double numerator_1= dot_product(wire.rc_firstder(phi),wire.rc_firstder(phi_p));
                double denominator_1=sqrt(r.distance(r_p)*r.distance(r_p)+wire.get_a()*wire.get_a()/sqrt(M_E));
                double integrand_1=numerator_1/denominator_1;

                //integral 2
                double numerator_2=dot_product(wire.rc_firstder(phi),wire.rc_firstder(phi));
                double denominator_2=sqrt(2*(1-cos(phi-phi_p))*numerator_2+wire.get_a()*wire.get_a()/sqrt(M_E));
                double integrand_2=numerator_2/denominator_2;

                //total integral
                return integrand_1-integrand_2;
            });
            double result_inner = gsl_integration_glfixed(inner, phi, phi+2*M_PI,
                                                          gsl_integration_glfixed_table_alloc(grid_size));
            return result_inner;
        });
        result_2D = mu_0/(4*M_PI)*gsl_integration_glfixed(outer, 0, 2*M_PI,
                                                      gsl_integration_glfixed_table_alloc(grid_size));

        //returning the final result
        return result_2D+result_1D;
    }
}