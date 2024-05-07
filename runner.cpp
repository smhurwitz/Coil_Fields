#include <iostream>
#include <functional>
#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include "src/l_calcs.cpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// REGULARIZED BIOT-SAVART LAW CONVERGENCE //////////////////////////////////////////////

/*
 * Plots the value of the regularized Biot-Savart Law across the length of the first HSX coil for several
 * different choices of quadrature grid points. The first column represents phi, the next columns represent the
 * value of the regularized Biot-Savart law at the correspond choice of quadrature grid points, and the final
 * column represents the converged solution.
 */
void biot_conv(double a, int N, string name, valarray<int> n_points_mod, valarray<int> n_points_unmod){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    Wire w = Wire::hsx(a, 1e6, 1);

    for(int i=0; i<N; i++){
        double phi = 2*i*M_PI/N;
        fprintf(fp, "%e ", phi);
        for (int n_pts : n_points_mod){
            double bx=b_reg(w, phi, 0, 2, n_pts);
            double by=b_reg(w, phi, 1, 2, n_pts);
            double bz=b_reg(w, phi, 2, 2, n_pts);
            double b = sqrt(bx*bx+by*by+bz*bz);
            fprintf(fp, "%e ", b);
        }

        for (int n_pts : n_points_unmod){
            double bx=b_reg(w, phi, 0, 1, n_pts);
            double by=b_reg(w, phi, 1, 1, n_pts);
            double bz=b_reg(w, phi, 2, 1, n_pts);
            double b = sqrt(bx*bx+by*by+bz*bz);
            fprintf(fp, "%e ", b);
        }

        double bx=b_reg(w, phi, 0, 3,-1);
        double by=b_reg(w, phi, 1, 3,-1);
        double bz=b_reg(w, phi, 2, 3,-1);
        double b = sqrt(bx*bx+by*by+bz*bz);
        fprintf(fp, "%e \n", b);
    }
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// CIRCULAR FIT COMPARISONS /////////////////////////////////////////////////////

/*
 * This method compares the circular fit method from Garren & Chen (column 3) to the circular fit method we obtained,
 * i.e., a choice of 0 quadrature points in the modified Biot-Savart law (column 2). The data is plotted on the first
 * HSX coil across its length. Column 1 represents phi. Column 4 represents the converged solution.
 */
void circ_fit_comp(double a, double I, string name, int N){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    Wire w = Wire::hsx(a, I, 1);

    for(int i=0; i<N; i++){
        double phi = 2*i*M_PI/N;
        fprintf(fp, "%e ", phi);
        double bx=b_reg(w, phi, 0, 2, 0);
        double by=b_reg(w, phi, 1, 2, 0);
        double bz=b_reg(w, phi, 2, 2, 0);
        double b = sqrt(bx*bx+by*by+bz*bz);
        fprintf(fp, "%e ", b);

        b=mu_0*I*w.kappa(phi)/(4*M_PI)*(log(8/(a*w.kappa(phi)))-0.75);
        fprintf(fp, "%e ", b);

        bx=b_reg(w, phi, 0, 3,-1);
        by=b_reg(w, phi, 1, 3,-1);
        bz=b_reg(w, phi, 2, 3,-1);
        b = sqrt(bx*bx+by*by+bz*bz);
        fprintf(fp, "%e \n", b);
    }
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAGNETIC FIELD CONTOURS //////////////////////////////////////////////////////
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
            cout << i*N+j << endl;
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

            Point p(s, theta, phi, coil);

            double bx = b_1D(p,0);
            double by = b_1D(p,1);
            double bz = b_1D(p,2);
            double modb = sqrt(bx * bx + by * by + bz * bz);

            fprintf(fp, "%e %e %e\n", x, y, modb);
        }
    }
    fclose(fp);
}

/*
 * This method outputs max|B| on a coil over phi∈[0,2\pi] with a choice of
 * N grid points. It is saved to a txt file of the name 'name.txt' in the format
 * [phi, max|B|_hifi, max|B|_semianalytic,max|B|_fullyanalytic], where phi and max|B| are both column vectors.
 */
void maxB_plots(Wire coil, int N, string name, double simplex=1e-3){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=138; i<N; i++){
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

/*
 * This method outputs the self-force of the first HSX coil against inverse aspect ratio at phi. The
 * aspect ratios are described by R, a_min, and a_max. The number of points to be plotted
 * is described by N.
 */
void force_vs_asp_ratio(int N, string name,double asp_rat_min, double asp_rat_max, double phi, int axis){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=0; i<N; i++){
        double asp_rat = asp_rat_min*pow(10,(log10(asp_rat_max)-log10(asp_rat_min))*i/(N-1));
        Wire w = Wire::hsx(0.326955182*asp_rat, 1, 1);
        double force = f_1D(w, phi, axis, 3);
        fprintf(fp, "%e %e \n", asp_rat, force);
    }
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// INDUCTANCE //////////////////////////////////////////////////////////////

/*
 * This method outputs the self-inductance of the first HSX coil against inverse aspect ratio. The
 * aspect ratios are described by R, a_min, and a_max. The number of points to be plotted
 * is described by N. The absolute and relative errors apply to the self-inductance calculations.
 * A good error is espabs=1e-5 and epsrel=1e-7.
 */
void L_plots(int N, string name, double epsabs, double epsrel, double asp_rat_min, double asp_rat_max){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=0; i<N; i++){
        double asp_rat = asp_rat_min*pow(10,(log10(asp_rat_max)-log10(asp_rat_min))*i/(N-1));
        Wire w = Wire::hsx(0.326955182*asp_rat, 1, 1);
        double inductance = L(w, 1e-8,1e-6);
        fprintf(fp, "%e %e \n", asp_rat, inductance);
    }
    fclose(fp);
}

/*
 * This method outputs the self-inductance of the first HSX coil against inverse aspect ratio. The
 * aspect ratios are described by R, a_min, and a_max. The number of points to be plotted
 * is described by N. The absolute and relative errors apply to the self-inductance calculations.
*/
void L_2D_plots(int N, string name, double epsabs, double epsrel, double asp_rat_min, double asp_rat_max){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    for(int i=0; i<N; i++){
        double asp_rat = asp_rat_min*pow(10,(log10(asp_rat_max)-log10(asp_rat_min))*i/(N-1));
        Wire w = Wire::hsx(0.326955182*asp_rat, 1, 1);
        double inductance = L_2D(w, epsrel, epsabs);
        fprintf(fp, "%e %e \n", asp_rat, inductance);
    }
    fclose(fp);
}

/*
 * This method exports to a .txt file the relative error in self-inductance versus number of Gauss-Legendre grid points. Column 1 is the
 * number of grid points per side, column 2 is the modified result, and column 3 is the unmodified result.
 */
void L_2D_convergence_plots(string name, int N_hi, Wire w){
    FILE *fp;
    name = name + ".txt";
    fp = fopen(name.c_str(),"w");
    double L_converged = L_2D_quadrature(w,400,2);
    for(int i=1; i<=N_hi; i++){
        double L_mod = L_2D_quadrature(w,i,2);
        double L_unmod = L_2D_quadrature(w,i,1);
        fprintf(fp, "%i %e %e \n", i, abs(L_mod-L_converged)/L_converged, abs(L_unmod-L_converged)/L_converged);
        cout << i << endl;
    }
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// TIMING COMPARISONS /////////////////////////////////////////////////////
/*
 * This method exports to a .txt file a set of (x,y,z) grid points on which to
 * evaluate the magnetic field
 */
void gen_grid(string name, int N, Wire w){
    FILE *fp;
    name = name + ".pts";
    fp = fopen(name.c_str(),"w");
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
                for(int k=0; k<N; k++){
                double u = w.get_a()*(2*i/N-1);
                double v = w.get_a()*(2*j/N-1);
                double s = sqrt(u*u + v*v);
                double theta = atan2(v, u);
                double phi = 2*k*M_PI/N;
                if (s <= w.get_a()){
                    Point p(s, theta, phi, w);
                    fprintf(fp, "%f\t%f\t%f\n", p.get_x(), p.get_y(), p.get_z());
                }
            }
        }
    }
    fclose(fp);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Gives some example code for computing data on important quantities such as self-force and self-inductance.
 */
int main(){
    // auto begin = std::chrono::high_resolution_clock::now(); //tracking how long code runs


    // // Calculations for an elliptical torus
    // int N = 1e0;
    // for (int n=0; n<N; n++){
    //     double a = 0.01;
    //     Wire hsx = Wire::hsx(a, 1e6, 1);
    //     Point p(a, 0, 0, hsx);
    //     // double bx = b_1D(p, 0);
    //     // double by = b_1D(p, 1);
    //     // double bz = b_1D(p, 2);
    //     double n_points = 1000;
    //     double bx = b_1D(p, 0, 2, n_points);
    //     double by = b_1D(p, 1, 2, n_points);
    //     double bz = b_1D(p, 2, 2, n_points);
    //     double modb = sqrt(bx * bx + by * by + bz * bz);
    //     cout << modb << endl;
    //     cout << std::to_string(p.get_x()) + ", "  + std::to_string(p.get_y()) + ", " +  std::to_string(p.get_z()) << endl;
    // }




    // double bx = b(p, 0, 1e-3, 1e-5);
    // double by = b(p, 1, 1e-3, 1e-5);
    // double bz = b(p, 2, 1e-3, 1e-5);

//    EXAMPLE CODE (note that length/(2*pi)=R=0.326955182 for the first HSX coil):

//    circ_fit_comp(0.326955182/100,1e6,"circ_fit",600); //comparison between circular fit methods

    gen_grid("grid", 20, Wire::hsx(0.01, 1e6, 1)); //generate grid for timing comparison

    //tracking how long code ran:
    // auto end = std::chrono::high_resolution_clock::now();
    // auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    // printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
    // printf("Time per field evaluation: %.4e seconds.\n", elapsed.count() * 1e-9 / N);
}