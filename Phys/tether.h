#ifndef TETHER_H
#define TETHER_H

//#include "force.h"
#include "q3d/edge.h"
QT_BEGIN_NAMESPACE
namespace Physics {
class Tether;
}
QT_END_NAMESPACE


class Physics::Tether
{
public:
    double _k,_scale;
    Tether(double k);
    void eval(Geometry::Edge *e);
    double calE(Geometry::Edge *e);
};
#endif // SPRINGFORCE_H
