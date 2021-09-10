#include "tether.h"

Physics::Tether::Tether(double k):
    _k(k),_scale(1)
{

}
void Physics::Tether::eval(Geometry::Edge *e) {
    double fx=0,fy=0,fz=0;
    double dx=DR(e->_vid1,e->_vid2,0);
    double dy=DR(e->_vid1,e->_vid2,1);
    double dz=DR(e->_vid1,e->_vid2,2);
    double dr=std::sqrt(dx*dx+dy*dy+dz*dz);
    double rl=e->_restLength;
    double cl0=rl*1.15*_scale;
    double cl1=rl*0.85/_scale;
    double max=rl*1.33*_scale;
    double min=rl*0.67/_scale;
    if(dr>max||dr<min){
        std::cout<<dr<<"\t"<<max<<"\t"<<min<<std::endl;
    }
    assert(dr<=max&&dr>=min);

    if(dr>cl0){
    //std::cout<<rl<<"\t"<<cl0<<"\t"<<dr<<std::endl;
        double dm=max-dr;
        double Idm=1/dm;
        double cl=cl0-dr;
        double Icl=1/cl;
        double constant=-_k*std::exp(Icl)*(1/(cl*cl*dm)+1/(dm*dm));
        fx=constant*dx/dr;
        fy=constant*dy/dr;
        fz=constant*dz/dr;
    }else
    if(dr<cl1){
        double dm=dr-min;
        double Idm=1/dm;
        double cl=dr-cl1;
        double Icl=1/cl;
        double constant=_k*std::exp(Icl)*(-1/(cl*cl*dm)+1/(dm*dm));
        fx=constant*dx/dr;
        fy=constant*dy/dr;
        fz=constant*dz/dr;
    }

    e->_vid1->_force->add(fx,fy,fz);
    e->_vid2->_force->add(-fx,-fy,-fz);

}
double Physics::Tether::calE(Geometry::Edge *e){
    double dx=DR(e->_vid1,e->_vid2,0);
    double dy=DR(e->_vid1,e->_vid2,1);
    double dz=DR(e->_vid1,e->_vid2,2);
    double dr=std::sqrt(dx*dx+dy*dy+dz*dz);
    double rl=e->_restLength;
    double cl0=rl*1.15*_scale;
    double cl1=rl*0.85/_scale;
    double max=rl*1.33*_scale;
    double min=rl*0.67/_scale;
    double ret=0;
    if(dr>cl0){
    //std::cout<<rl<<"\t"<<cl0<<"\t"<<dr<<std::endl;
        double dm=max-dr;
        double cl=cl0-dr;
        double Icl=1/cl;
        ret=_k*std::exp(Icl)/dm;

    }else
    if(dr<cl1){
        double dm=dr-min;
        double cl=dr-cl1;
        double Icl=1/cl;
        ret=_k*std::exp(Icl)/dm;

    }
    return ret;
}
