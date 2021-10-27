#ifndef QUAD_H
#define QUAD_H

#include <Phys/bead.h>
#include <Phys/springforce.h>
QT_BEGIN_NAMESPACE
namespace Physics {
class Quad;
}
QT_END_NAMESPACE


using Physics::VecD3d;
class Physics::Quad:public QObject
{
private:
    Physics::VecD3d *_temp;

public:
    double _lens[4];
    Physics::Bead **_beads,*_mb;
    Physics::SpringForce *_sf;
    Quad(Physics::Bead*,double);
    void addNewCandidate(Physics::Bead*);
    void eval();

};

#endif // QUAD_H
