#ifndef BENDINGPARAMETERS_H
#define BENDINGPARAMETERS_H
#include <qobject.h>
#include "utils/classes.h"
#include "q3d/beadinfo.h"
#include "Phys/vecd3d.h"
#include "Phys/tensor2.h"
#include "q3d/edge.h"

#include "Phys/bendingenergy.h"


#define fixIndex(j,list) ((j)>=0?((j)%list->size()):list->size()+(j))
#define getTheOtherBead(e,b) (e->_vid1==b?e->_vid2:e->_vid1)

static const int LeviC[3][3][3] = {{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};
QT_BEGIN_NAMESPACE
namespace Physics {
class BendingParameters;
}

QT_END_NAMESPACE

class Physics::BendingParameters
{


public:
    double _T;
    double _triA;
    double chiM,chiP;
    double _cot1,_cot2;
    double _l,_lsq,_phi;
    Physics::VecD3d *_dxP,*_normal;
    Geometry::BeadInfo *_bi;
    Physics::Tensor2 *_tensor;
    Geometry::Triangle * _tri;



    BendingParameters(Geometry::BeadInfo *bi);

    void update(Geometry::BeadInfo*);
    void updateBendingParameters();
};

#endif // BENDINGPARAMETERS_H
