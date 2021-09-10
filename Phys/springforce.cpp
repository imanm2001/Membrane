#include "springforce.h"

Physics::SpringForce::SpringForce(double k):_k(k)
{

}
double Physics::SpringForce::eval(Geometry::Edge *e) {
    double dx=DR(e->_vid1,e->_vid2,0);
    double dy=DR(e->_vid1,e->_vid2,1);
    double dz=DR(e->_vid1,e->_vid2,2);
    double dr=std::sqrt(dx*dx+dy*dy+dz*dz);
    //double k=_k*std::exp(4*std::fabs(e->_restLength-dr)/e->_restLength);
    double k=_k;

    dx*=k*(e->_restLength-dr)/dr;
    dy*=k*(e->_restLength-dr)/dr;
    dz*=k*(e->_restLength-dr)/dr;
    e->_vid1->_force->add(dx,dy,dz);
    e->_vid2->_force->add(-dx,-dy,-dz);
    double a=dr/e->_restLength;
    dr=e->_restLength-dr;

    return _k*dr*dr/2;

}
