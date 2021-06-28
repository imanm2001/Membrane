#ifndef SPRINGFORCE_H
#define SPRINGFORCE_H

//#include "force.h"
#include "q3d/edge.h"

QT_BEGIN_NAMESPACE
namespace Physics {
class SpringForce2;
}
QT_END_NAMESPACE


class Physics::SpringForce2
{
public:
    double _k,_len;
    SpringForce2(double k,double len);
    void eval(Physics::Bead *b1,Physics::Bead *b2);
};
#endif // SPRINGFORCE_H
