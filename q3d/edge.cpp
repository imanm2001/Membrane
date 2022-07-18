#include "edge.h"

void Geometry::Edge::location(Physics::VecD3d* ret){
    ret->setValues(_vid1->_coords);
    ret->add(_vid2->_coords);
    ret->multConst(0.5);
}
Geometry::Edge::Edge(QObject * p,Physics::BeadInfo* vid1,Physics::BeadInfo* vid2):_restLengthScale(1)
{
    _vid1=vid1;
    _vid2=vid2;
    assert(_vid1!=vid2);
    _tris[0]=nullptr;
    _tris[1]=nullptr;
    this->setParent(p);
    //FIXME

    if(dynamic_cast<Geometry::Mesh*>(p) != nullptr){
        Geometry::Mesh *u=((Geometry::Mesh*)p);
        u->addNewEdge(this);
    }
    //updateRL();

    double dx=DR(_vid1,_vid2,0);
    double dy=DR(_vid1,_vid2,1);
    double dz=DR(_vid1,_vid2,2);
    _restLength=std::sqrt(dx*dx+dy*dy+dz*dz);
    _vid1->_connections->append(this);
    _vid2->_connections->append(this);
    //_restLength=1/19.0;


}
inline void Geometry::Edge::updateRL(){

    double dx=DR(_vid1,_vid2,0);
    double dy=DR(_vid1,_vid2,1);
    double dz=DR(_vid1,_vid2,2);
    //_restLength=std::sqrt(dx*dx+dy*dy+dz*dz);
}
bool Geometry::Edge::equal(Physics::BeadInfo* b1,Physics::BeadInfo * b2){
    return (_vid1==b1&&_vid2==b2)||(_vid1==b2&&_vid2==b1);
}
bool Geometry::Edge::equal(Geometry::Edge * e){
    return (_vid1==e->_vid1&&_vid2==e->_vid2)||(_vid1==e->_vid2&&_vid2==e->_vid1);
}
int Geometry::Edge::triIndex(Geometry::Triangle* t){
    int ret=-1;
    if(_tris[0]==t){
        ret=0;
    }else if(_tris[1]==t){
        ret=1;
    }
    return  ret;
}
Physics::BeadInfo* Geometry::Edge::commonVertex(Geometry::Edge *e2){
    Physics::BeadInfo* ret=nullptr;
    if(_vid1==e2->_vid1||_vid1==e2->_vid2){
        ret=_vid1;
    }else  if(_vid2==e2->_vid1||_vid2==e2->_vid2){
        ret=_vid2;
    }
    return ret;
}
void Geometry::Edge::print(){
    std::cout<<_vid1<<"\t"<<_vid2<<std::endl;
}

bool Geometry::Edge::contain(Physics::BeadInfo* vid){
    return _vid1==vid||_vid2==vid;
}
Geometry::Edge::~Edge(){

}
