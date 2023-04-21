#include <valarray>
#include <stdexcept>
#include <math.h>
using namespace std;

/* Returns the vector norm of a vector. */
double vector_norm(valarray<double> v){
    double sum = 0;
    for(double x : v){sum += pow(x,2);}
    return sqrt(sum);
}

/* Performs a cross product on the 3-vectors u and v. */
valarray<double> cross_product(valarray<double> u, valarray<double> v){
    if(u.size() != 3 || v.size() != 3){throw std::invalid_argument( "vector is incorrect size");} //make sure vectors are length 3.
    valarray<double> uxv={0,0,0};
    uxv[0]=u[1]*v[2]-u[2]*v[1];
    uxv[1]=u[2]*v[0]-u[0]*v[2];
    uxv[2]=u[0]*v[1]-u[1]*v[0];
    return uxv;
}

/* Performs a dot product on the 3-vectors u and v. */
double dot_product(valarray<double> u, valarray<double> v){
    if(u.size() != 3 || v.size() != 3){throw std::invalid_argument( "vector is incorrect size");} //make sure vectors are length 3.
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}
