#include "bead.h"

Physics::Bead::Bead(QObject* p,Physics::VecD3d* coord,double GAMMA,int id):QObject(p),ID(id),_coords(coord),_GAMMA(GAMMA){
    _force=new VecD3d();
    _force->zero();
}
Physics::Bead::Bead(Bead *parent):Bead(parent->parent(),new VecD3d(parent->_coords),parent->_GAMMA,parent->ID){


}
void Physics::Bead::update(double dt){
    //GAMMA * D=KbT
    _coords->_coords[0]+=_force->_coords[0]*dt/_GAMMA;
    _coords->_coords[1]+=_force->_coords[1]*dt/_GAMMA;
    _coords->_coords[2]+=_force->_coords[2]*dt/_GAMMA;
    _force->zero();
}
