#ifndef SPRINGFORCE_H
#define SPRINGFORCE_H

//#include "force.h"
#include "q3d/edge.h"

QT_BEGIN_NAMESPACE
namespace Physics {
class SpringForce;
}
QT_END_NAMESPACE


class Physics::SpringForce
{
public:
    double _k;
    SpringForce(double k);
    double eval(Geometry::Edge *e,double *dist);
    double calE(Geometry::Edge *e);
};
#endif // SPRINGFORCE_H
