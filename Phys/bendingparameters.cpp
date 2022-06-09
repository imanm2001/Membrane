
#include "bendingparameters.h"

Physics::BendingParameters::BendingParameters(Geometry::BeadInfo *bi):_bi(bi)
{

    _dxP=new VecD3d();
    _tensor=new Tensor2();
    //_normal=new VecD3d(0,1,0);
    _normal=nullptr;
    _lsq=_T=_l=0;
    _phi=chiM=chiP=0;
    _tri=nullptr;

}

void Physics::BendingParameters::update(Physics::BeadInfo* bj){
    assert(_bi!=bj);

    _dxP->setValues(_bi->_coords);
    _dxP->sub(bj->_coords);
    _lsq=_dxP->dot(_dxP);
    _l=std::sqrt(_lsq);
    if(_l<=0||_l!=_l){

        _bi->_coords->print();
        bj->_coords->print();
        std::cout<<_l<<std::endl;
    }
    /*
    assert(_l>0);
    assert(_l==_l);
    assert(!isinf(_l));*/
    IERROR(_l>0);
    IERROR(_l==_l);
    IERROR(!isinf(_l));


}

