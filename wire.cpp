#include <iostream>
#include <valarray>
#include <functional>
#include <stdexcept>
#include <math.h>
#include "vectors.cpp"
#include "integration.cpp"
using namespace std;

/*
 * This class stores information on a particular wire loop described by the set of Fourier coefficients
 * for each Cartesian component of r_c(phi), i.e r_c^(x)=Sum_n(A_n*cos(n*phi+phi_n)). This class also provides an
 * analytical representation at any point of rc(phi), its first and second derivatives, its curvature,
 * and Frenet-Serret unit vectors.
 */
class Wire {
private:
    double I; //current
    double a; //minor radius
    int order; //number of elements (including n=0 term) in the Fourier series
    valarray<double> xs_coeff, xc_coeff, ys_coeff, yc_coeff, zs_coeff, zc_coeff; //fourier coefficients

public:
    Wire(){};
    Wire(double I, double a, int order, valarray<double> xc_coeff, valarray<double> xs_coeff,
                                        valarray<double> yc_coeff, valarray<double> ys_coeff,
                                        valarray<double> zc_coeff, valarray<double> zs_coeff) {
        this->I = I;
        this->a = a;
        this->order = order;

        // Check to make sure that all the coefficient vectors are of the expected order. If not, throw an error.
        if(xc_coeff.size() != order || yc_coeff.size() != order || zc_coeff.size() != order ||
           xs_coeff.size() != order || ys_coeff.size() != order || zs_coeff.size() != order){
            throw std::invalid_argument( "vectors are incorrect size");
        }

        this->xc_coeff=xc_coeff; this->xs_coeff=xs_coeff;
        this->yc_coeff=yc_coeff; this->ys_coeff=ys_coeff;
        this->zc_coeff=zc_coeff; this->zs_coeff=zs_coeff;
    }

    /*
     * This function returns the Wire representation of a torus of radii R and a and current I.
     */
    static Wire torus(double R, double a, double I){
        return Wire(I, a, 2, {0,R}, {0,0}, {0,0},
                                   {0,R}, {0,0}, {0,0});
    }

    /*
     * This function returns the Wire representation of a HSX coil obtained from SIMSOPT. There are 6 coils provided,
     * and coil_number=1,...,6 can select which coil to use. The SIMSOPT representation is further detailed in
     * https://simsopt.readthedocs.io/en/latest/simsopt_user.geo.html#simsopt.geo.CurveXYZFourier.
     */
    static Wire hsx(double a, double I, int coil_number){
        valarray<double> simsopt_coefficients;
        if(coil_number == 1){
            simsopt_coefficients = {1.45009861e+00, -2.47414872e-01, -4.46196554e-02,  8.63915549e-03,
                                    -1.69234274e-02, -1.78265183e-02, -1.27439499e-02, -1.71361226e-03,
                                    -1.32923421e-03, -1.00582328e-04, -3.36107311e-03, -2.77955139e-03,
                                    7.96480116e-04,  6.26375236e-04,  1.29940775e-04, -7.09183134e-04,
                                    -2.88738126e-04, -3.98533232e-04,  7.78052075e-06, -7.81615347e-05,
                                    -2.00984745e-04, -5.14353694e-04, -2.65380474e-04, -2.24969687e-04,
                                    7.80576092e-05, -2.06331897e-04, -1.24938525e-04, -3.02233154e-04,
                                    6.68072202e-05, -7.66618316e-05,  1.13622400e-04, -1.29823846e-04,
                                    3.90715806e-05,  8.13648141e-02, -1.67423912e-02, -1.78507840e-01,
                                    -1.29080809e-02, -2.79157508e-03, -5.10964909e-03,  1.17200057e-02,
                                    1.08618043e-03,  3.85166510e-03, -2.69721424e-03,  7.70845580e-03,
                                    2.05279160e-03, -7.69885459e-06, -1.63125193e-03,  2.78503206e-03,
                                    2.12289016e-03,  1.17337799e-03,  2.33820331e-04, -2.42942408e-04,
                                    8.06624087e-04,  8.25543637e-04,  8.43239011e-04, -6.39529164e-04,
                                    8.18676546e-06,  3.34609368e-06,  4.45393873e-04, -1.93725100e-04,
                                    -7.75940044e-05, -2.60611778e-04,  3.06908552e-07,  2.96890093e-05,
                                    -2.34814636e-05, -8.23928877e-05,  5.21477002e-02, -7.88660646e-02,
                                    3.08534912e-01, -2.97820617e-03, -6.34432891e-04, -1.14097111e-02,
                                    1.73079693e-02,  4.64977901e-03,  5.58022381e-03, -7.68821305e-04,
                                    1.19348515e-03,  1.18988153e-03,  1.06888912e-03, -2.84674973e-04,
                                    8.69405639e-04,  2.96494415e-04,  7.15670205e-04, -5.05642926e-05,
                                    4.56171080e-04, -7.80768755e-05,  4.16485375e-04, -1.88586814e-05,
                                    1.60017831e-04, -1.72992119e-04,  2.04604496e-04, -4.80898483e-07,
                                    9.77190383e-05, -5.86528582e-05,  1.90233949e-05, -1.02867931e-05,
                                    8.66851714e-06,  3.72016159e-05, -6.15318270e-05};
        }
        else if (coil_number == 2){
            simsopt_coefficients = {1.36424714e+00, -1.97423333e-01, -6.66973535e-02,  1.52791255e-02,
                                    -2.17500802e-02, -1.48160502e-02, -2.52097382e-02, -5.57515715e-03,
                                    -5.53378881e-03, -2.59078419e-03, -1.25840243e-02, -2.50069826e-03,
                                    2.10702848e-03, -2.15894351e-03, -7.87210194e-04, -1.19037871e-03,
                                    3.40412328e-04, -8.14908668e-04,  1.26291058e-03, -5.52339381e-04,
                                    -2.48193560e-05, -1.07944430e-04,  1.27718234e-04, -1.38147484e-04,
                                    3.53788909e-04, -4.31615785e-05, -2.32501263e-04, -2.84969651e-05,
                                    1.89391558e-04, -5.43715475e-05,  1.28350062e-05, -4.75493871e-05,
                                    2.04701195e-05,  2.29188189e-01, -8.91359845e-02, -1.78615403e-01,
                                    -1.72919938e-02,  4.37221068e-03, -1.78382026e-03,  1.55155529e-02,
                                    3.52333192e-03,  4.94881954e-03,  5.92282762e-03,  1.38579387e-02,
                                    5.92354443e-04, -1.74587089e-03,  1.88468853e-03,  2.24742973e-03,
                                    1.47016006e-03, -7.89025608e-04, -9.01892787e-04, -8.23524207e-04,
                                    1.06818258e-03, -1.79651019e-04, -6.67804687e-04, -2.03180188e-04,
                                    1.19239020e-04,  5.76480245e-05,  5.15762333e-05,  1.55306333e-04,
                                    -9.97763414e-06,  1.70430002e-04,  1.01859013e-04,  1.65652555e-05,
                                    9.40603478e-05,  1.02878530e-04,  1.41555575e-01, -1.29258934e-01,
                                    3.19995324e-01, -5.51796083e-03, -7.82300500e-03, -1.70029736e-02,
                                    1.55638284e-02,  4.14293534e-03, -9.64996736e-04, -1.35328034e-03,
                                    2.86285888e-03, -1.21504507e-04, -1.28806034e-03, -9.89263454e-04,
                                    1.25458485e-03, -2.14810947e-04,  8.62355574e-04, -1.16672498e-03,
                                    7.22112342e-04,  5.83017758e-04,  1.03744965e-03, -4.65764067e-04,
                                    4.75236470e-04,  2.48557751e-04,  3.63496583e-04,  1.74671655e-04,
                                    9.85916511e-05, -2.85794267e-05,  3.61991748e-05,  1.33687296e-04,
                                    -1.85898219e-04,  1.07891854e-05, -5.93430717e-05};
        }
        else if(coil_number == 3){
            simsopt_coefficients={1.22708494e+00, -1.25651192e-01, -1.16143188e-01,  2.75615911e-02,
                                  -2.37812160e-02, -1.96548747e-02, -3.53380310e-02, -3.65681970e-03,
                                  -6.58550835e-03, -1.22511077e-02, -1.33261568e-02, -3.09894355e-04,
                                  -1.24154975e-03, -1.07180301e-03,  1.29163706e-03, -1.08823524e-03,
                                  1.65589387e-03,  5.38492180e-04, -3.87768025e-04, -3.48118117e-05,
                                  1.22470076e-03, -8.97943902e-04,  1.03089599e-04, -1.80739749e-04,
                                  -7.78705862e-04, -1.45311288e-04,  5.01045666e-04, -1.61784768e-04,
                                  -1.41031815e-04, -8.15404470e-05, -1.69856324e-04,  1.84347676e-04,
                                  1.72192121e-04,  3.58693196e-01, -1.30160611e-01, -1.91713627e-01,
                                  -2.45611317e-02,  8.42360384e-03,  7.42700197e-03,  1.52386896e-02,
                                  -4.12227727e-03,  1.33308025e-03,  1.68305580e-02,  9.76035993e-03,
                                  -2.80188734e-03,  8.84399682e-04,  1.81418962e-03,  1.44031990e-03,
                                  3.92182273e-04, -2.37196825e-03, -1.20634073e-03,  1.90957982e-03,
                                  -5.96413374e-04, -1.20196865e-03,  7.74119125e-04,  1.05051855e-04,
                                  9.46197330e-05,  3.05596642e-05,  1.48889531e-04, -5.72143213e-05,
                                  2.56834332e-04, -3.08733306e-04,  2.62408300e-05,  2.97710227e-04,
                                  -1.72699288e-04,  6.60675834e-06,  1.81213313e-01, -1.59565412e-01,
                                  3.14525445e-01, -1.11703225e-02, -1.38116927e-02, -1.87871481e-02,
                                  7.27964954e-03, -5.97653252e-03, -9.68080119e-03, -8.20241116e-04,
                                  -7.37261378e-04, -2.82497228e-03, -2.96417482e-03, -1.53683101e-03,
                                  3.16248413e-03,  5.35727670e-04, -5.21292431e-05, -1.05935004e-03,
                                  3.03032758e-03,  4.58482543e-04,  8.27732038e-04,  3.01811910e-04,
                                  4.69768267e-04,  1.81446325e-04,  3.13666959e-04, -2.08752404e-04,
                                  1.49774525e-05,  2.11000101e-04, -3.79091029e-04, -1.38490375e-04,
                                  4.92073596e-05,  8.92236814e-06, -7.93877358e-05};
        }
        else if(coil_number == 4){
            simsopt_coefficients={1.08292784e+00, -1.02050612e-01, -1.53569440e-01,  2.92444765e-02,
                                  -4.20433041e-02, -4.56409163e-02, -1.15786671e-02, -3.56629771e-03,
                                  -1.30267875e-02, -5.89741736e-03,  8.35825041e-03,  1.03231046e-03,
                                  3.29642602e-03, -1.28663860e-03,  1.36242233e-03,  1.61208658e-03,
                                  -1.96891085e-03,  9.30973468e-04, -1.40282521e-03, -1.63513991e-04,
                                  9.57052477e-04, -3.84917296e-04,  3.23539280e-04,  6.58276326e-05,
                                  -6.78898186e-04, -4.81781248e-06, -2.03422384e-04,  3.21774508e-04,
                                  9.18995317e-05, -2.43900370e-04,  2.70663826e-04, -1.07258942e-04,
                                  2.93956652e-05,  4.75889920e-01, -1.90426756e-01, -1.73330323e-01,
                                  -3.24747593e-02,  2.53626030e-02,  1.87789530e-02, -9.48783682e-03,
                                  3.11142777e-03,  1.71501953e-02,  4.20949244e-03, -1.23586191e-02,
                                  6.72082539e-04,  1.30776618e-04,  1.23920239e-03, -4.11116693e-03,
                                  1.68628200e-03,  2.08967935e-03, -1.87904272e-03, -5.74563237e-04,
                                  7.31380620e-04,  7.45148927e-05, -4.46093132e-05, -8.78707504e-04,
                                  1.81181099e-05,  3.62065231e-04, -1.80272028e-04,  5.83054557e-04,
                                  -4.10688614e-05, -1.15880205e-04,  7.36136624e-05, -8.51103040e-05,
                                  1.54198969e-04, -5.39442917e-05,  1.55508566e-01, -1.06148430e-01,
                                  2.98899498e-01, -2.12012381e-02, -1.44391859e-02, -1.79890942e-02,
                                  6.13455825e-03, -1.70929528e-02,  5.31332401e-05, -4.39109349e-03,
                                  2.61074658e-03, -4.05032891e-03,  5.20008142e-03,  1.13792620e-03,
                                  1.44564100e-03,  3.51860676e-03, -5.07885114e-05, -3.22920115e-04,
                                  -1.21437327e-03,  1.08834478e-03,  2.19536703e-05, -9.32514460e-04,
                                  -1.25462954e-03, -1.43093972e-05,  1.73267067e-04, -2.25552990e-04,
                                  4.47067424e-04,  2.07565370e-05,  1.35973896e-04, -2.41390134e-04,
                                  9.87493164e-05,  2.82490157e-04, -1.06009162e-04};
        }
        else if(coil_number == 5){
            simsopt_coefficients={9.35577932e-01, -3.57387618e-02, -2.08341527e-01,  3.42863624e-02,
                                  -2.11405934e-02, -5.66641489e-02, -9.45030671e-03, -8.35005446e-03,
                                  -1.18974258e-02, -3.67845887e-03, -1.44986983e-03,  7.49641613e-03,
                                  -1.45826456e-03,  3.38318754e-03,  3.02439848e-03, -2.50963357e-03,
                                  6.59622802e-04, -1.33908080e-03,  8.86955719e-04, -7.95855306e-04,
                                  -8.20640481e-04,  2.88003693e-05, -9.42642078e-04,  8.23889101e-05,
                                  3.85553751e-04,  4.19892836e-04,  1.28087250e-05, -2.00943782e-04,
                                  9.40777105e-05, -2.00335639e-04,  3.62006231e-04,  9.61640428e-05,
                                  2.95309526e-04,  5.67322186e-01, -1.72176022e-01, -2.01685673e-01,
                                  -3.45034203e-02,  2.57269182e-02,  2.69999277e-02, -1.16151220e-02,
                                  1.01040607e-02,  1.92636582e-02,  5.80701809e-03, -4.24265749e-03,
                                  4.42380068e-04,  1.65763644e-03, -9.60086618e-04, -2.70922285e-03,
                                  1.22469913e-03,  1.31957638e-03,  3.47188878e-04, -5.65271548e-04,
                                  4.11041838e-05, -8.06798291e-04, -4.16072346e-04, -1.65188182e-04,
                                  1.83912059e-04, -4.31219727e-04,  9.34854065e-05,  1.46415242e-04,
                                  -1.14892875e-04,  1.23081553e-04,  1.63309610e-04,  8.43214378e-05,
                                  5.81867849e-05, -1.51186461e-04,  1.05897435e-01, -1.62956524e-01,
                                  1.89333521e-01, -6.53013323e-03, -3.41768394e-02, -3.62088491e-02,
                                  -2.49338982e-03, -1.62472470e-02, -5.18110393e-03, -5.63889605e-03,
                                  1.69276159e-03,  2.60952856e-03,  5.24121382e-03,  8.88110505e-04,
                                  9.98329154e-04,  2.07288284e-03,  1.62286991e-03,  6.65381552e-04,
                                  -2.07084499e-03,  2.90915452e-04, -5.81647148e-05, -1.05714656e-03,
                                  8.05616530e-04,  8.18285464e-05,  6.41109009e-04, -4.14265295e-04,
                                  7.24708636e-05,  7.25051252e-05, -5.61224070e-05,  2.21328984e-04,
                                  1.55975466e-04,  1.09722857e-04, -1.55480154e-04};
        }
        else if(coil_number == 6){
            simsopt_coefficients={7.82360135e-01,  8.87865022e-02, -2.09682904e-01,  2.25491599e-02,
                                  2.56785115e-02, -4.78858950e-02, -2.40623194e-02, -1.04433524e-02,
                                  -7.82706885e-03, -4.94701207e-05, -9.50318507e-03,  6.05639498e-03,
                                  -1.94048836e-03,  1.30832259e-03,  7.68221447e-04, -1.25232687e-03,
                                  -3.20574050e-04, -1.64216232e-03, -1.11618397e-04, -1.17993418e-03,
                                  -3.11834479e-04, -5.53316670e-04,  4.51926811e-04,  3.01708005e-04,
                                  1.58813215e-05,  2.71874178e-04, -4.93474192e-05, -2.20564558e-04,
                                  2.29888930e-04, -2.34467228e-04,  3.52454955e-04, -4.13566021e-05,
                                  -4.42748113e-05,  6.60448993e-01, -1.06066394e-01, -2.29973452e-01,
                                  -3.60747291e-02,  2.58303671e-02,  3.76516696e-02,  8.90159824e-03,
                                  8.20441478e-03,  1.58436391e-02,  1.04261762e-02,  2.78345476e-04,
                                  2.21692618e-03, -2.80132299e-03,  1.79846602e-03, -8.73467278e-04,
                                  -7.62568932e-04,  1.41029433e-04, -1.47785960e-03, -6.40028314e-04,
                                  -6.28526867e-04, -3.46390407e-04, -1.86017369e-04, -3.43212260e-04,
                                  3.97273540e-04,  2.00648512e-04,  2.29942796e-04,  1.33108380e-04,
                                  3.22412092e-04, -1.00673337e-04,  3.30945863e-04, -1.54792372e-04,
                                  1.95968513e-04, -1.19505375e-04,  3.90080533e-02, -2.02770466e-01,
                                  3.50058435e-02,  1.42962453e-02, -1.64470132e-02, -5.57677669e-02,
                                  -2.26738119e-02, -9.64549991e-03, -1.14870805e-02, -2.09091988e-03,
                                  2.13780807e-04, -4.64330065e-04,  5.75018214e-03,  2.22179083e-03,
                                  3.09563979e-03,  5.69582904e-04,  3.73967806e-04,  2.18941978e-04,
                                  -8.16675486e-04, -5.42788003e-04, -2.74092136e-04, -2.49302577e-04,
                                  3.05609776e-04,  2.26438768e-04,  7.06503581e-05, -8.33218939e-06,
                                  1.81235851e-04, -4.50224910e-04,  5.13681100e-04, -3.30205845e-04,
                                  3.36808319e-04,  1.02784668e-04, -2.26230572e-05};
        }
        else {
            {throw std::invalid_argument( "coil_number must be an integer between 1 and 6");}
        }
        int order = 1 + (simsopt_coefficients.size() / 3 - 1) / 2; //order of the fourier series
        valarray<double> xs(order), xc(order), ys(order), yc(order), zs(order), zc(order);
        for (int i = 0; i < order; i++){
            xc[i] = simsopt_coefficients[2*i];
            yc[i] = simsopt_coefficients[2 * order - 1 + 2 * i];
            zc[i] = simsopt_coefficients[(2 * order - 1) * 2 + 2 * i];
            if(i>0) {
                xs[i] = simsopt_coefficients[2 * i - 1];
                ys[i] = simsopt_coefficients[2 * order - 2 + 2 * i];
                zs[i] = simsopt_coefficients[(2 * order - 1) * 2 + 2 * i - 1];
            }
            else{
                xs[i]=0;
                ys[i]=0;
                zs[i]=0;
            }
        }
        return Wire(I, a, order, xc, xs, yc, ys, zc, zs);
    }

    double get_a(){return a;}
    double get_I(){return I;}

    /*
     * This function gives the value of r_c(phi) in terms of the Cartesian components, [r_c^(x), r_c^(y), rc^(z)].
     */
    valarray<double> rc(double phi){
        valarray<double> rc = {0, 0, 0};
        for(int i=0; i < order; i++){
            rc[0]+= xc_coeff[i] * cos(i * phi) + xs_coeff[i] * sin(i * phi);
            rc[1]+= yc_coeff[i] * cos(i * phi) + ys_coeff[i] * sin(i * phi);
            rc[2]+= zc_coeff[i] * cos(i * phi) + zs_coeff[i] * sin(i * phi);
        }
        return rc;
    }

    /*
     * This function gives the value of the first derivative of rc(phi) in terms of the Cartesian components,
     * [r'_c^(x), r'_c^(y), r'_c^(z)].
     */
    valarray<double> rc_firstder(double phi){
        valarray<double> r_c_der = {0,0,0};
        for(int i=0; i < order; i++){
            r_c_der[0]+=i*(-xc_coeff[i]*sin(i*phi)+xs_coeff[i]*cos(i*phi));
            r_c_der[1]+=i*(-yc_coeff[i]*sin(i*phi)+ys_coeff[i]*cos(i*phi));
            r_c_der[2]+=i*(-zc_coeff[i]*sin(i*phi)+zs_coeff[i]*cos(i*phi));
        }
        return r_c_der;
    }

    /* This function returns the value of the vector norm of the first dervative of rc(phi), |r'_c(phi)|. */
    double rc_norm_firstder(double phi){
        return vector_norm(rc_firstder(phi));
    }

    /*
     * This function gives the value of the second derivative of rc(phi) in terms of the Cartesian components,
     * [r''_c^(x), r''_c^(y), r''_c^(z)].
     */
    valarray<double> r_c_second_der(double phi){
        valarray<double> r_c_der = {0,0,0};
        for(int i=0; i < order; i++){
            r_c_der[0]+=-i*i*(xc_coeff[i]*cos(i*phi)+xs_coeff[i]*sin(i*phi));
            r_c_der[1]+=-i*i*(yc_coeff[i]*cos(i*phi)+ys_coeff[i]*sin(i*phi));
            r_c_der[2]+=-i*i*(zc_coeff[i]*cos(i*phi)+zs_coeff[i]*sin(i*phi));
        }
        return r_c_der;
    }

    /*
     * This function returns the value of the first derivative of the vector norm of the first dervative of rc(phi),
     * d(|r'_c(phi)|)/dphi.
     */
    double rc_der_normfirstder(double phi){
        double sum1 =0, sum2=0, sum3=0, sum4=0, sum5=0, sum6=0;
        for(int i=0; i < order; i++){
            sum1+=i*(-xc_coeff[i]*sin(i*phi)+xs_coeff[i]*cos(i*phi));
            sum2+=i*i*(xc_coeff[i]*cos(i*phi)+xs_coeff[i]*sin(i*phi));
            sum3+=i*(-yc_coeff[i]*sin(i*phi)+ys_coeff[i]*cos(i*phi));
            sum4+=i*i*(yc_coeff[i]*cos(i*phi)+ys_coeff[i]*sin(i*phi));
            sum5+=i*(-zc_coeff[i]*sin(i*phi)+zs_coeff[i]*cos(i*phi));
            sum6+=i*i*(zc_coeff[i]*cos(i*phi)+zs_coeff[i]*sin(i*phi));
        }
        return (sum1*sum2+sum3*sum4+sum5*sum6) / rc_norm_firstder(phi);
    }

    /*
     * This function returns the value of the vector norm of the curvature.
     */
    double kappa(double phi){
        double norm = rc_norm_firstder(phi);
        valarray<double> first_der = rc_firstder(phi);
        valarray<double> second_der = r_c_second_der(phi);
        return vector_norm(cross_product(first_der,second_der))/pow(norm,3);
    }

    /*
     * This function gives the value of the first unit vector, e1. This vector points along r'(phi). The Cartesian
     * components are returned.
     */
    valarray<double> e1(double phi) {
        valarray<double> der_r_c = rc_firstder(phi);
        return der_r_c / vector_norm(der_r_c);
    }

    /*
     * This function gives the value of the second unit vector, e2. This vector points along the curvature. The
     * Cartesian components are returned.
     */
    valarray<double> e2(double phi) {
        valarray<double> second_der = r_c_second_der(phi);
        double norm_second_der=vector_norm(second_der);
        return second_der/norm_second_der;
    }

    /*
     * This function gives the value of the third unit vector, e3. This vector is orthogonal to
     * e1 and e2. The Cartesian components are returned.
     */
    valarray<double> e3(double phi) {return cross_product(e1(phi), e2(phi));}

    /*
     * This function returns the length of the wire by integrating the quantity |dr_c/dphi| over phi.
     */
    double length(){
        auto len_integrand = make_gsl_function([&](double phi){
            return rc_norm_firstder(phi);
        });

        double result, abserr;
        gsl_integration_qags(len_integrand, 0, 2*M_PI, 1e-8, 1e-8,
                             100,IntegrationWorkspace(100),&result, &abserr);
        return result;
    }
};