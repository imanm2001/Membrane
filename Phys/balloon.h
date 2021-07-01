#ifndef BALLOON_H
#define BALLOON_H
#include <QObject>
#include "utils/classes.h"
#include "q3d/wavefrontobj.h"
#include "Phys/springforce.h"
#include "math.h"
#include "surfacewithphysics.h"
QT_BEGIN_NAMESPACE
namespace Physics {
class Balloon;
}
QT_END_NAMESPACE

class Physics::Balloon:public SurfaceWithPhysics
{
private:
    Geometry::WaveFrontObj * _sphere;
    double _dt;
    Physics::VecD3d *_temp;
    Physics::SpringForce *_sf;
    int _step;
    QRandomGenerator* _rand;
public:
    double _l,_lsq,_totalA;

    Balloon(double dt);
    Qt3DRender::QGeometryRenderer * mesh();
    void update();
    void updateVIN();
    double calE();
    double calE(int);
};

#endif // BALLOON_H
