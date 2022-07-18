#include "springforce.h"

Physics::SpringForce::SpringForce(double k):_k(k)
{

}
double Physics::SpringForce::eval(Geometry::Edge *e,double* dist) {
    double dx=DR(e->_vid1,e->_vid2,0);
    double dy=DR(e->_vid1,e->_vid2,1);
    double dz=DR(e->_vid1,e->_vid2,2);
    double dr=std::sqrt(dx*dx+dy*dy+dz*dz);
    dist[0]=dr;
    //double k=_k*std::exp(4*std::fabs(e->_restLength-dr)/e->_restLength);
    double k=_k;
    double RL=e->_restLength*e->_restLengthScale;
    double DR=RL-dr;
    dx*=k*(DR)/dr;
    dy*=k*(DR)/dr;
    dz*=k*(DR)/dr;
    e->_vid1->_force->add(dx,dy,dz);
    e->_vid2->_force->add(-dx,-dy,-dz);
    double a=dr/e->_restLength;


    return _k*DR*DR/2;

}
double Physics::SpringForce::calE(Geometry::Edge *e) {
    double dx=DR(e->_vid1,e->_vid2,0);
    double dy=DR(e->_vid1,e->_vid2,1);
    double dz=DR(e->_vid1,e->_vid2,2);
    double dr=std::sqrt(dx*dx+dy*dy+dz*dz);
    dr=e->_restLength*e->_restLengthScale-dr;
    return _k*dr*dr/2;

}
