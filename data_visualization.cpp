#include <iostream>
#include <functional>
#include <stdexcept>
#include <stdio.h>
#include <math.h>

#include "f_calcs.cpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// MAGNETIC FIELD ///////////////////////////////////////////////////////////
/*
 * This high-fidelity method outputs |B| across a slice phi. The data is sent to a text file of the name 'name.txt'
 * in the form [x, y, |B|], where x, y, and |B| are column vectors. The data is printed on an NxN grid.
 */
void contours(Wire coil, double phi, int N, string name){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=0; i<N; i++) {
        for (int j = 0; j < N; j++) {
            double a = coil.get_a();
            a *= 1.2; // giving a bit of extra room on the plot
            double x = -a + 2 * a * i / (N - 1);
            double y = -a + 2 * a * j / (N - 1);
            double s = sqrt(x * x + y * y);
            double theta = atan2(y, x);
            Point p(s, theta, phi, coil);
            double bx = b(p, 0, 1e-3, 1e-5);
            double by = b(p, 1, 1e-3, 1e-5);
            double bz = b(p, 2, 1e-3, 1e-5);
            double modb = sqrt(bx * bx + by * by + bz * bz);
            fprintf(fp, "%e %e %e\n", x, y, modb);
        }
    }
}

/*
 * This low-fidelity method outputs |B| across a slice phi. The data is sent to a text file of the name 'name.txt'
 * in the form [x, y, |B|], where x, y, and |B| are column vectors. The data is printed on an NxN grid.
 */
void contours_1D(Wire coil, double phi, int N, string name){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=0; i<N; i++) {
        for (int j = 0; j < N; j++) {
            double a = coil.get_a();
            double x = -a*1.2 + 2 * a*1.2 * i / (N - 1);
            double y = -a*1.2 + 2 * a*1.2 * j / (N - 1);
            double s = sqrt(x * x + y * y);
            double theta = atan2(y, x);

            /*
             * We don't have analytical formulas for b_extra outisde the body of the wire. Purely for plotting purposes
             * we will try to hide contours outside of the wire by setting the force to NaN.
             */
            if (s>a){fprintf(fp, "%e %e %e\n", x, y, nan(""));}
            else{
                Point p(s, theta, phi, coil);

                double bx = b_1D(p,0);
                double by = b_1D(p,1);
                double bz = b_1D(p,2);
                double modb = sqrt(bx * bx + by * by + bz * bz);

                fprintf(fp, "%e %e %e\n", x, y, modb);
            }
        }
    }
    fclose(fp);
}

/*
 * This method outputs max|B| on a coil over phi∈[0,2\pi] along an axis (0=x, 1=y, 2=z) with a choice of
 * N grid points. It is saved to a txt file of the name 'name.txt' in the format
 * [phi, max|B|_hifi, max|B|_semianalytic,max|B|_fullyanalytic], where phi and max|B| are both column vectors.
 */
void maxB_plots(Wire coil, int N, string name, double simplex=1e-3){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=0; i<N; i++){
        double phi = 2*i*M_PI/N;
        double bhi=max_modB(coil,phi,simplex);
        double b_semianalytic = max_modB_semianalytic(coil,phi, 1e-6);
        double b_fullyanalytic = max_modB_fullyanalytic(coil,phi);
        fprintf(fp, "%e %e %e %e\n", phi, bhi, b_semianalytic, b_fullyanalytic);
    }
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// FORCE ////////////////////////////////////////////////////////////////
/*
 *  This method outputs the force per unit length on a coil over phi∈[0,2\pi] along an axis {x,y,z}
 *  corresponding to an int {0,1,2}, respectively. It is saved to a txt file of the name 'name.txt' in the format
 *  [phi, hifi, lofi], where phi, hifi, and lofi are column vectors. The number of grid points is specified by 'N'. The
 *  choice of the low-fidelity method is specific by 'key', where 0 corresponds to the 1D integral method and 1
 *  corresponds to the circular approximation method.
 */
void force_plots(Wire coil, int axis, int N, string name, int key, double epsrel, double epsabs){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=0; i<N; i++){
        double phi = 2*i*M_PI/N;
        double fhi = f(coil, phi, axis, epsrel, epsabs);
        double flo;
        if(key==0){flo = f_1D(coil, phi, axis, 3, -1);}
        else{flo = f_circular(coil,phi,axis);}
        fprintf(fp, "%e %e %e\n", phi, fhi, flo);
    }
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(){
    auto begin = std::chrono::high_resolution_clock::now();
//    EXAMPLES:

//    force_plots(Wire::hsx(0.00326955182, 1e6, 1), 0, 200, "f", 0,0.1*1e-12, 0.1*1e-12);
//    contours(Wire::hsx(0.00326955182, 1e6, 1), 0, 60, "exact01");
    maxB_plots(Wire::hsx(0.00326955182, 1e6, 1),1, "maxB", 1e-3);
//    Wire w = Wire::hsx(0.00326955182/100, 1e6, 1);
//cout << b(Point(0.00326955182/100,0,0,w),1,1e-4,1e-4) << endl;
//cout << b_1D(Point(0.00326955182/100,0,0,w),1);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);

}