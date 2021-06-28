#include "vecd3d.h"

Physics::VecD3d::VecD3d():VecD3d(0,0,0)
{

}
void Physics::VecD3d::subVec(VecD3d * v2,VecD3d * res){
    res->_coords[0]=_coords[0]-v2->_coords[0];
    res->_coords[1]=_coords[1]-v2->_coords[1];
    res->_coords[2]=_coords[2]-v2->_coords[2];
}
Physics::VecD3d::VecD3d(double x,double y,double z){
    _coords[0]=x;
    _coords[1]=y;
    _coords[2]=z;
}
Physics::VecD3d::VecD3d(double vals[3]){
    _coords[0]=vals[0];
    _coords[1]=vals[1];
    _coords[2]=vals[2];
}
Physics::VecD3d::VecD3d(QVector3D* v){
    _coords[0]=v->x();
    _coords[1]=v->y();
    _coords[2]=v->z();
}
Physics::VecD3d::VecD3d(VecD3d* v){
    setValues(v);
}
void Physics::VecD3d::setValues(VecD3d* v){
    _coords[0]=v->_coords[0];
    _coords[1]=v->_coords[1];
    _coords[2]=v->_coords[2];
}
double Physics::VecD3d::len(){
    return std::sqrt(_coords[0]*_coords[0]+_coords[1]*_coords[1]+_coords[2]*_coords[2]);
}
double Physics::VecD3d::dot(VecD3d*v){
    return  _coords[0]*v->_coords[0]+_coords[1]*v->_coords[1]+_coords[2]*v->_coords[2];
}
void Physics::VecD3d::cross(VecD3d* a,VecD3d* r){
    r->_coords[0]=_coords[1]*a->_coords[2]-_coords[2]*a->_coords[1];
    r->_coords[1]=_coords[2]*a->_coords[0]-_coords[0]*a->_coords[2];
    r->_coords[2]=_coords[0]*a->_coords[1]-_coords[1]*a->_coords[0];
}
void Physics::VecD3d::print(){
    std::cout<<_coords[0]<<"\t"<<_coords[1]<<"\t"<<_coords[2]<<"\r\n";
}
void Physics::VecD3d::zero(){
    _coords[0]=_coords[1]=_coords[2]=0;
}
void Physics::VecD3d::add(double dx,double dy,double dz){
_coords[0]+=dx;
_coords[1]+=dy;
_coords[2]+=dz;
}
void Physics::VecD3d::sub(double dx,double dy,double dz){
_coords[0]-=dx;
_coords[1]-=dy;
_coords[2]-=dz;
}
void Physics::VecD3d::sub(Physics::VecD3d *v){
    _coords[0]-=v->_coords[0];
    _coords[1]-=v->_coords[1];
    _coords[2]-=v->_coords[2];
}

void Physics::VecD3d::nomilize(){
    double l=len();
    _coords[0]/=l;
    _coords[1]/=l;
    _coords[2]/=l;

}
void Physics::VecD3d::add(VecD3d v){
    _coords[0]+=v._coords[0];
    _coords[1]+=v._coords[1];
    _coords[2]+=v._coords[2];

}
void Physics::VecD3d::multConst(double s){
    _coords[0]*=s;
    _coords[1]*=s;
    _coords[2]*=s;

}
Physics::VecD3d::~VecD3d(){

}
void Physics::VecD3d::setValues(double x,double y,double z){
    _coords[0]=x;
    _coords[1]=y;
    _coords[2]=z;
}
void Physics::VecD3d::debug(){
    assert(_coords[0]==_coords[0]);
    assert(_coords[1]==_coords[1]);
    assert(_coords[2]==_coords[2]);
}
