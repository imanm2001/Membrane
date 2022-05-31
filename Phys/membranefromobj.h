#ifndef MEMBRANEFROMOBJ_H
#define MEMBRANEFROMOBJ_H
#include <QObject>

#include "utils/classes.h"
#include "q3d/wavefrontobj.h"
#include "Phys/springforce.h"
#include "math.h"
#include "surfacewithphysics.h"
#include "Phys/tether.h"
#include "Phys/CTS.h"
#include <gsl/gsl_linalg.h>
QT_BEGIN_NAMESPACE
namespace Physics {
class MembraneFromObj;
}
QT_END_NAMESPACE
class Physics::MembraneFromObj:public SurfaceWithPhysics
{
protected:
    void capture(double *data,int len);
private:
    const int _NUMSAVEDATA=3;
    double _saveData[3];
    VecD3d  **_CTSvs1,**_CTSvs2;
    Geometry::BeadInfo **_CTSbeads;
    gsl_matrix *_strainMatrix;
    gsl_vector *_strainDirection,*_strainDirection2;
    double _initalArea,_pBE,_appliedF,_py,_cfy,_APV,_pA,_dtF,_odtF,_maxR;
    double _pE[8];
    Physics::CTS *_cts;
    Geometry::WaveFrontObj * _disc;
    double _dt,_tension,_ptension,_ptension2,_TIE,_radialForce,_MRad,_FSign;
    Physics::VecD3d *_temp,*_temp2,*_temp3;
    Physics::SpringForce *_sf;
    Physics::Tether *_tet;
    int _step,_Rind,_pstep;
    Physics::BeadInfo* _cb;
    QVector<Geometry::BeadInfo*> *_border,*_xprofile,*_updatable;
    double _RFA,_APFA;
    void updateBeads(QVector<Geometry::BeadInfo*> *,double P,double kappa,double dt);
    QRandomGenerator* _rand;
    double _frad;
    double _kappaFactor,_radiusFactor;
    void testStrain();
    void findThefourthVertex(Geometry::Triangle*,Geometry::Triangle*,Geometry::BeadInfo**);
    Geometry::Triangle* _tempTri;
public:
    MembraneFromObj(double dt);
    Qt3DRender::QGeometryRenderer * mesh();
    void update();
    void updateVIN();
    double calE();
    double calE(int);
    double calStrain();
    double calStrain2D();
    double calStrain2D2();
    void calRect(VecD3d**,double*,double*);
    void Project2D(VecD3d**,VecD3d*,VecD3d*,VecD3d*);



};

#endif // MEMBRANEFROMOBJ_H
