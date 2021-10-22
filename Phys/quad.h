#ifndef QUAD_H
#define QUAD_H

#include <Phys/bead.h>
class Quad
{
public:
    double _l1,_l2,_l3,_l4;
    Physics::Bead *_b1,*_b2,*_b3,*_b4;
    Quad(Physics::Bead*);

};

#endif // QUAD_H
