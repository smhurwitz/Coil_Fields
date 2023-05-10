

#include "b_calcs.cpp"
#include <functional>
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// HI-FI ///////////////////////////////////////////////////////////
struct ForceData{double phi; Wire wire; int axis;};

/*
 * This is the integrand for the function 'f'.
 */
static int f_integrand(const int *ndim, const cubareal xx[],
                       const int *ncomp, cubareal ff[], void *data){
    // Note that the Cuhre integration routine integrates over a unit cube, so each integration coordinate (below)
    // must be multiplied by a certain factor to ensure the integral has the proper domain.
    ForceData *u = (struct ForceData*)data;
    int axis = u->axis;
    Wire wire = u->wire;
    double phi = u->phi;
    double s = xx[0]*wire.get_a();
    double theta = xx[1]*2*M_PI;
    double s_p= xx[2] * wire.get_a();
    double theta_p= xx[3] * 2 * M_PI+theta;
    double phi_p= xx[4] * 2 * M_PI+phi;

    Point r(s, theta, phi, wire);
    Point r_p(s_p, theta_p, phi_p, wire);

    double denominator = pow(r_p.distance(r)+1e-100, 3);
    double num1 = s * s_p * (1 - wire.kappa(phi) * s * cos(theta)) * (1 - wire.kappa(phi_p) * s_p * cos(theta_p));
    valarray<double> e1 = wire.e1(phi);
    valarray<double> numerator = cross_product(e1,cross_product(wire.rc_firstder(phi_p), r.vec_distance(r_p)));

    ff[0] = (wire.get_a()*wire.get_a())*(num1 * numerator[axis] / denominator);
}

/*
 *  This high-fidelity method returns the value of dF/dl the magnetic field at a slice phi and along an axis {x,y,z}
 *  corresponding to an int {0,1,2}, respectively.
 */
double f(Wire wire, double phi, int axis, double epsrel, double epsabs){
    int nregions, neval, fail;
    cubareal integral[1], error[1], prob[1];
    struct ForceData data;
    data.phi = phi;
    data.axis = axis;
    data.wire = wire;

    Cuhre(5, 1, f_integrand, &data, 1, epsrel, epsabs, 0,0,
          1000000000, 0,nullptr, nullptr,&nregions, &neval, &fail, integral,
          error, prob);
    cout << neval << endl;
    double prefactors =(2*mu_0*pow(wire.get_I(),2)/(wire.get_a()*wire.get_a()))/(wire.get_a()*wire.get_a());
    return prefactors*integral[0];
}

/*
 * This is the integrand for the function 'f_modified'.
 */
static int f_modified_integrand(const int *ndim, const cubareal xx[],
                                const int *ncomp, cubareal ff[], void *data){
    ForceData *u = (struct ForceData*)data;
    int axis = u->axis;
    Wire wire = u->wire;
    double phi = u->phi;
    double s = xx[0]*wire.get_a();
    double theta = xx[1]*2*M_PI;
    double s_p= xx[2] * wire.get_a();
    double theta_p= xx[3] * 2 * M_PI+theta;
    double phi_p= xx[4] * 2 * M_PI+phi;

    Point r(s, theta, phi, wire);
    Point r_p(s_p, theta_p, phi_p, wire);

    double denominator = pow(r_p.distance(r) + 1e-100, 3);
    double num1 = s * s_p * (1 - wire.kappa(phi) * s * cos(theta)) *
                  (1 - wire.kappa(phi_p) * s_p * cos(theta_p))
                  * wire.rc_norm_firstder(phi_p);
    valarray<double> e1 = wire.e1(phi);
    valarray<double> e2 = wire.e2(phi);
    valarray<double> e3 = wire.e3(phi);
    valarray<double> e1_p = wire.e1(phi_p);
    valarray<double> numerator = cross_product(e1,cross_product(e1_p, {r.dx(r_p), r.dy(r_p), r.dz(r_p)}));

    valarray<double> e1_pf = sin(phi_p - phi) * e2 + cos(phi_p - phi) * e1;
    valarray<double> e2_pf = cos(phi_p - phi) * e2 - sin(phi_p - phi) * e1;
    valarray<double> e3_pf = e3;
    valarray<double> rc_pf = wire.rc(phi)+1/wire.kappa(phi)*((1-cos(phi_p - phi))*e2 +sin(phi_p - phi)*e1);

    valarray<double> r_pf = rc_pf + s_p * cos(theta_p) * e2_pf + s_p * sin(theta_p) * e3_pf;
    valarray<double> r_vec = r.get_vector_form();
    double denominator2 = pow(vector_norm(r_vec - r_pf) + 1e-100, 3);
    double num12 = s * s_p * (1 - wire.kappa(phi) * s * cos(theta)) *
                   (1 - wire.kappa(phi) * s_p * cos(theta_p)) / (wire.kappa(phi));
    valarray<double> numerator2 = cross_product(e1, cross_product(e1_pf, r_vec - r_pf));

    ff[0] = (wire.get_a()*wire.get_a())*(num1 * numerator[axis] / denominator - num12 * numerator2[axis] / denominator2);
    return nan("");
}

/*
 *  This semi-high-fidelity method returns the value of dF/dl the magnetic field at a slice phi and along an axis {x,y,z}
 *  corresponding to an int {0,1,2}, respectively. The modified hi-fi method refers to the act of subtracting a known
 *  function off from the integrand of the original hi-fi method in order to make it easier to numerically integrate.
 */
double f_modified(Wire wire, double phi, int axis, double epsrel, double epsabs){
    int nregions, neval, fail;
    cubareal integral[1], error[1], prob[1];
    struct ForceData data;
    data.phi = phi;
    data.axis = axis;
    data.wire = wire;

    Cuhre(5, 1, f_modified_integrand, &data, 1, epsrel, epsabs, 0,0,
          1000000000, 0,nullptr, nullptr,&nregions, &neval, &fail, integral,
          error, prob);

    valarray<double> analytic_result = -mu_0*pow(wire.get_I(),2)*wire.kappa(phi)/(4*M_PI)*(-0.75+
            log(8/(wire.kappa(phi) * wire.get_a()))) * wire.e2(phi); //analytic force for a torus
    double prefactors =(2*mu_0*pow(wire.get_I(),2)/(wire.get_a()*wire.get_a()))/(wire.get_a()*wire.get_a());
    return prefactors*integral[0]+analytic_result[axis];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// LO-FI ///////////////////////////////////////////////////////////
/*
 * This low-fidelity method returns the value of dF/dl the magnetic field at a slice phi and along an axis {x,y,z}
 * corresponding to an int {0,1,2}, respectively. This method evaluates the 1D integral in one of three ways:
 *      (1) applying Gauss-Legendre integration to the unmodified integral,
 *      (2) applying Gauss-Legendre integration to the modified integral, or
 *      (3) using an adaptive routine to accurately evaluate the integral.
 * The method of choice can be selected with the integer 'key', where (1=unmodified, 2=modified, 3=adaptive). If the user
 * selects a Gauss-Legendre routine, they may specify the number of grid points 'n_points'. Otherwise, this variable
 * is unnecessary.
 */
double f_1D(Wire wire, double phi, int axis, int key, int n_points=50){
    valarray<double> vector_current = wire.get_I()* wire.e1(phi);
    if(axis==0){
        return vector_current[1] * b_reg(wire, phi, 2, key, n_points) -
               vector_current[2] * b_reg(wire, phi, 1, key, n_points);
    }
    else if(axis==1){
        return vector_current[0] * b_reg(wire, phi, 2, key, n_points) -
               vector_current[2] * b_reg(wire, phi, 0, key, n_points);
    }
    else{
        return vector_current[0] * b_reg(wire, phi, 1, key, n_points) -
               vector_current[1] * b_reg(wire, phi, 0, key, n_points);
    }
}

/*
 * This low-fidelity circular approximation method returns the value of dF/dl the magnetic field at a slice phi and
 * along an axis {x,y,z} corresponding to an int {0,1,2}, respectively
 */
double f_circular(Wire wire, double phi, int axis){
    valarray<double> analytic_result = -mu_0*pow(wire.get_I(),2)*wire.kappa(phi)/(4*M_PI)*(-0.75+
            log(8/(wire.kappa(phi) * wire.get_a()))) * wire.e2(phi);
    return analytic_result[axis];
}
