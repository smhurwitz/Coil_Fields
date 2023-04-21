#include <valarray>
#include <functional>
#include <math.h>
#include "wire.cpp"
using namespace std;

/* This Class describes a point in 3D space based on the Frenet-Serret coordinate system (s,theta,phi; Wire). */
class Point{
private:
    double s, theta, phi; //Frenet-Serret representation of the point
    double x, y, z; //alternative Cartesian representation of the point
    valarray<double> vec_rep; //Cartesian vector representation, i.e. [x, y, z]
    Wire wire;

public:
    Point(){};
    Point(double s, double theta, double phi, Wire wire){
        this->s=s;
        this->theta=theta;
        this->phi=phi;
        this->wire = wire;

        valarray<double> rc = wire.rc(phi);
        valarray<double> e2 = wire.e2(phi);
        valarray<double> e3 = wire.e3(phi);

        vec_rep = rc + s*cos(theta)*e2 + s*sin(theta)*e3;
        x=vec_rep[0]; y=vec_rep[1]; z=vec_rep[2];
    }

    /*
     * Return the Cartesian distance between this Point and the Point 'x'.
     */
    double distance(Point x){
        valarray<double> rc = wire.rc(phi); valarray<double> rc_x = wire.rc(x.get_phi());
        valarray<double> e2 = wire.e2(phi); valarray<double> e2_x = wire.e2(x.get_phi());
        valarray<double> e3 = wire.e3(phi); valarray<double> e3_x = wire.e3(x.get_phi());

        valarray<double> dr = rc-rc_x;
        valarray<double> de2 = s*cos(theta)*e2-x.get_s()*cos(x.get_theta())*e2_x;
        valarray<double> de3 = s*sin(theta)*e3-x.get_s()*sin(x.get_theta())*e3_x;

        return vector_norm(dr+de2+de3);
    }

    Wire get_wire(){return wire;}

    double get_s() const {return s;}
    double get_theta() const {return theta;}
    double get_phi() const {return phi;}

    valarray<double> get_vector_form(){return vec_rep;}
    double get_x(){return x;}
    double get_y(){return y;}
    double get_z(){return z;}
    double dx(Point p){return x-p.get_x();}
    double dy(Point p){return y-p.get_y();}
    double dz(Point p){return z-p.get_z();}
};