#include "triangle.h"

Geometry::Triangle::Triangle(QObject* p,int id,Geometry::Edge *e1,Geometry::Edge *e2,Geometry::Edge *e3,double c)
    :QObject(p),_ID(id),_curvater(c)
{
    setParent(p);
    _location=new Physics::VecD3d();
    _e[0]=e1;
    _e[1]=e2;
    _e[2]=e3;
    addToEdge(e1);
    addToEdge(e2);
    addToEdge(e3);


    auto v1=e1->commonVertex(e3);
    auto v2=e1->commonVertex(e2);
    auto v3=e3->commonVertex(e2);
#ifdef QT_DEBUG
    assert(v1&&v2&&v3&&v1!=v2&&v1!=v3&&v2!=v3);
#endif
    _v[0]=v1;
    _v[1]=v2;
    _v[2]=v3;
    _norm=new VecD3d();
    getNormal();
    /*
    if(getNormal()._coords[2]<0){
        auto vt=v2;
        v2=v1;
        v1=vt;
    }
    _v[0]=v1;
    _v[1]=v2;
    _v[2]=v3;*/


    //updateVIndecies();
}
Geometry::Triangle::Triangle(QObject *p, int id,Physics::BeadInfo* i0,Physics::BeadInfo* i1,Physics::BeadInfo* i2,double c):
    Triangle(p,id,new Geometry::Edge(p,i0,i1),new Geometry::Edge(p,i1,i2),new Geometry::Edge(p,i2,i0),c)
{

}
void Geometry::Triangle::addToEdge(Geometry::Edge *e){
#ifdef QT_DEBUG
    assert(e->_tris[1]==nullptr);
#endif
    if(e->_tris[0]==nullptr){
        e->_tris[0]=this;
    }else {
        e->_tris[1]=this;
    }
}
bool Geometry::Triangle::isCommonEdge(Geometry::Edge *e){
    return e==_e[0]||e==_e[1]||e==_e[2];
}
Geometry::Edge* Geometry::Triangle::commonEdge(Triangle *t){
    Geometry::Edge* ret=nullptr;
    if(t!=nullptr){
        for(int i=0;i<3;i++){
            if(isCommonEdge(t->_e[i])){
                ret=t->_e[i];
                break;
            }
        }
    }
    return ret;
}
/*
void Geometry::Triangle::updateVIndecies(){
    auto v0=_e[0]->_vid1;
    auto v1=_e[0]->_vid2;
    Physics::Bead* comV=_e[0]->commonVertex(_e[1]);
    auto v2=comV==_e[1]->_vid1?_e[1]->_vid2:_e[1]->_vid1;
    Physics::Bead * vs[3],*nb=nullptr;;
    vs[0]=v0;
    vs[1]=v1;
    vs[2]=v2;
    int oldIndex[2],index=0,i2=0;
    for(int i=0;i<3;i++){
        int k;
        if((k=includeVertex(vs[i]))>-1){
            oldIndex[index++]=k;

        }else{
            nb=vs[i];
            i2=i;

        }
    }
    for(int i=0;i<3;i++){
        if(i!=oldIndex[0]&&i!=oldIndex[1]){
            index=i;
            break;
        }
    }
    _v[index]=nb;
    getNormal();
#ifdef QT_DEBUG
    if(_norm._coords[2]<0){
        std::cout<<"Normal ERRROORRRR------"<<std::endl;
    }
    for(int i=0;i<3;i++){
        assert(_e[i]->_tris[0]!=_e[i]->_tris[1]);
    }
#endif
}*/

void Geometry::Triangle::updateVIndecies(){
    Physics::BeadInfo *newV=nullptr;
    Physics::BeadInfo * commV[2];
    int oldInd=-1;
    int cmvi=0;
    for(int i=0;i<3;i++){
        if((_e[0]->contain(_v[i])||_e[1]->contain(_v[i])||_e[2]->contain(_v[i]))){
            commV[cmvi++]=_v[i];
        }else{
#ifdef QT_DEBUG
            assert(oldInd==-1);
#endif
            oldInd=i;
        }
    }
#ifdef QT_DEBUG
    assert(oldInd!=-1);
    assert(cmvi==2);
#endif
    for(int i=0;i<3;i++){
        if(_e[i]->_vid1!=commV[0]&&_e[i]->_vid1!=commV[1]){
            newV=_e[i]->_vid1;
            break;
        }
        if(_e[i]->_vid2!=commV[0]&&_e[i]->_vid2!=commV[1]){
            newV=_e[i]->_vid2;
            break;
        }
    }
    /*
    if(oldInd==1){

        _v[1]=newV;
        newV=_v[0];
        _v[0]=_v[2];
        _v[2]=newV;

    }*/
    _v[oldInd]=newV;
#ifdef QT_DEBUG
    auto v1=_e[0]->commonVertex(_e[2]);
    auto v2=_e[0]->commonVertex(_e[1]);
    auto v3=_e[2]->commonVertex(_e[1]);

    assert(v1&&v2&&v3&&v1!=v2&&v1!=v3&&v2!=v3);

    for(int i=0;i<3;i++){
        if(_v[i]==nullptr){
            std::cout<<"NULLL"<<i<<"_"<<oldInd<<std::endl;
            assert(_v[i]!=nullptr);
        }
    }
#endif
    //getNormal();
#ifdef QT_DEBUG
    assert(newV!=nullptr);
/*
    if(_norm._coords[2]<0){
        std::cout<<"Normal ERRROORRRR------\t"<<oldInd<<std::endl;
    }*/
    for(int i=0;i<3;i++){
        assert(_e[i]->_tris[0]!=_e[i]->_tris[1]);
    }
#endif

}
int Geometry::Triangle::includeVertex(Physics::BeadInfo* b){

    for(int i=0;i<3;i++){
        if(_v[i]==b){
            return i;
        }
    }
    return -1;
}
Geometry::Edge ** Geometry::Triangle::edgeWithVertex(Physics::BeadInfo* vid){
    int index=0;
    Geometry::Edge **ret=(Geometry::Edge **)malloc(sizeof (Geometry::Edge*)*2);
    for(int i=0;i<3&&index<2;i++){
        if(_e[i]->contain(vid)){
            ret[index++]=_e[i];
        }
    }
    return ret;
}
Geometry::Edge * Geometry::Triangle::edgeWithVertexExclude(Physics::BeadInfo* vid,Geometry::Edge* ex){
    Geometry::Edge *ret=nullptr;
    assert(_e);
    for(int i=0;i<3;i++){
        assert(_e[i]);
        if(_e[i]!=ex&&_e[i]->contain(vid)){
            ret=_e[i];
            break;
        }
    }
    return ret;
}
Physics::VecD3d* Geometry::Triangle::getNormal(){
    double dx1=DR(_v[1],_v[0],0);
    double dy1=DR(_v[1],_v[0],1);
    double dz1=DR(_v[1],_v[0],2);

    double dx2=DR(_v[2],_v[0],0);
    double dy2=DR(_v[2],_v[0],1);
    double dz2=DR(_v[2],_v[0],2);

    double dx3=DR(_v[2],_v[1],0);
    double dy3=DR(_v[2],_v[1],1);
    double dz3=DR(_v[2],_v[1],2);

    _norm->_coords[0]=det(dy1,dy2,dz1,dz2);
    _norm->_coords[1]=det(dz1,dz2,dx1,dx2);
    _norm->_coords[2]=det(dx1,dx2,dy1,dy2);
    _area=_norm->len()/2;

    _norm->nomilize();

    #ifdef QT_DEBUG
    /*
    if(_norm._coords[2]<0){
        std::cout<<"NERROR\t"<<((Plane*)parent())->_tris.indexOf(this) <<std::endl;


    }*/
    #endif

    double dl1=std::sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
    double dl2=std::sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
    double dl3=std::sqrt(dx3*dx3+dy3*dy3+dz3*dz3);

    double dot1=dx1*dx2+dy1*dy2+dz1*dz2;
    double dot2=-(dx1*dx3+dy1*dy3+dz1*dz3);

    _angles[0]=std::acos(dot1/(dl1*dl2));
    _angles[1]=std::acos(dot2/(dl1*dl3));
    _angles[2]=PI-_angles[0]-_angles[1];
    _isObtuse=_angles[0]>PI2||_angles[1]>PI2||_angles[2]>PI2;

    return _norm;
}
void Geometry::Triangle::printEdges(){
    std::cout<<_e[0]<<"\t"<<_e[1]<<"\t"<<_e[2]<<std::endl;
}
double Geometry::Triangle::area(){
    return _area;
}
Physics::VecD3d* Geometry::Triangle::getLocation(){
    _location->setValues(_v[0]->_coords);
    _location->add(_v[1]->_coords);
    _location->add(_v[2]->_coords);
    _location->multConst(1/3.0);
    return _location;
}
int Geometry::Triangle::getVertexIndex(BeadInfo* bi){
    int ret=-1;
    if(bi==_v[0]){
        ret=0;
    }else if(bi==_v[1]){
        ret=1;
    }else if(bi==_v[2]){
        ret=2;
    }
    return ret;
}
bool Geometry::Triangle::contains(BeadInfo* bi){
    return _v[0]==bi||_v[1]==bi||_v[2]==bi;
}
