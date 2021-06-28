#ifndef MEMBRANEFROMOBJ_H
#define MEMBRANEFROMOBJ_H
#include <QObject>

#include "utils/classes.h"
#include "q3d/wavefrontobj.h"
#include "Phys/springforce.h"
#include "math.h"
#include "surfacewithphysics.h"
#include "Phys/tether.h"
QT_BEGIN_NAMESPACE
namespace Physics {
class MembraneFromObj;
}
QT_END_NAMESPACE
class Physics::MembraneFromObj:public SurfaceWithPhysics
{
protected:
    void capture();
private:
    double _initalArea,_pBE;
    Geometry::WaveFrontObj * _disc;
    double _dt,_tension;
    Physics::VecD3d *_temp;
    Physics::SpringForce *_sf;
    Physics::Tether *_tet;
    int _step;
    Physics::BeadInfo* _cb;
    QVector<Geometry::BeadInfo*> *_border,*_xprofile;
    void updateBeads(QVector<Geometry::BeadInfo*> *,double P,double kappa,double dt);
    QRandomGenerator* _rand;

public:
    MembraneFromObj(double dt);
    Qt3DRender::QGeometryRenderer * mesh();
    void update();
    void updateVIN();
    double calE();
    double calE(int);



};

#endif // MEMBRANEFROMOBJ_H
