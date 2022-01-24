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
    double _initalArea,_pBE,_appliedF,_py,_cfy,_APV,_pA;
    double _pE[8];

    Geometry::WaveFrontObj * _disc;
    double _dt,_tension,_ptension,_ptension2,_TIE,_radialForce,_MRad,_FSign;
    Physics::VecD3d *_temp,*_temp2;
    Physics::SpringForce *_sf;
    Physics::Tether *_tet;
    int _step,_Rind;
    Physics::BeadInfo* _cb;
    QVector<Geometry::BeadInfo*> *_border,*_xprofile,*_updatable;
    void updateBeads(QVector<Geometry::BeadInfo*> *,double P,double kappa,double dt);
    QRandomGenerator* _rand;
    double _frad;
    double _kappaFactor,_radiusFactor;

public:
    MembraneFromObj(double dt);
    Qt3DRender::QGeometryRenderer * mesh();
    void update();
    void updateVIN();
    double calE();
    double calE(int);



};

#endif // MEMBRANEFROMOBJ_H
