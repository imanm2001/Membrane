#include "springforce2.h"

Physics::SpringForce2::SpringForce2(double k,double len):_k(k),_len(len)
{

}
void Physics::SpringForce2::eval(Physics::Bead *b1,Physics::Bead *b2) {
    VecD3d *res=new VecD3d();
    b1->_coords->subVec(b2->_coords,res);
//    double dx=DR(b1,b2,0);
//    double dy=DR(b1,b2,1);
//    double dz=DR(b1,b2,2);
//    double dr=std::sqrt(dx*dx+dy*dy+dz*dz);
//    dx*=_k*(_len-dr)/dr;
//    dy*=_k*(_len-dr)/dr;
//    dz*=_k*(_len-dr)/dr;
    double dr=res->len();
    res->multConst(_k*(_len-dr)/dr);
    b1->_force->add(res);
    b2->_force->sub(res);

}
